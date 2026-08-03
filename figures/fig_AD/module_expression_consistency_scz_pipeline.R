# Of the modules that are expressed differently in control vs disease, what module genes are driving this?
# DEseq2 documentation: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library(DESeq2)
library(qs)
library(data.table)
library(AnnotationHub)
library(tidyverse)
library(cowplot)
options(bitmapType = 'cairo')

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/module_expression_consistency_fxns.R"))


##########################
# Run pipeline for SCZ (bulk megaset)
##########################
# (bulk megaset modules onto all brainSCOPE data)
# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/brainSCOPE/")
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/scz")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       sn_summary_path = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/",
                       return_obj = T,
                       plot = T,
                       comp_names = c("SCZ", "Con"),
                       header_names = c("brainSCOPE"))

# bulk megaset onto CMC
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/brainSCOPE_CMC/")
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_CMC/conVsSCZ/")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       sn_summary_path = copa_dir1,
                       return_obj = T,
                       plot = T,
                       comp_names = c("SCZ", "Con"),
                       header_names = c("brainSCOPE_CMC"))

# bulk megaset onto SZBD
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/brainSCOPE_SZBDMultiseq/")
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_SZBDMultiseq/conVsSCZ/")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       sn_summary_path = copa_dir1,
                       return_obj = T,
                       plot = T,
                       comp_names = c("SCZ", "Con"),
                       header_names = c("brainSCOPE_SZBD"))

##########################
# Make list of overlapping modules between cohorts (brainSCOPE cohorts)
##########################
paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_CMC"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_SZBDMultiseq"))
pathnames <- c("CMC", "SZBD")
 
# Load list of significant modules by subclass
countdf_unfilt <- lapply(paths, \(d){
  compname <- list.files(d)
  modlist <- lapply(compname, \(x){
    path <- file.path(d, x, "modlist_filtered.qs")
    if(file.exists(path)){
      df <- qread(file.path(d, x, "modlist_filtered.qs"))
      outlist <- mapply(\(a, b){
        if(nrow(a) > 0){
          out <- a[, 1:2]
          out$type = b
          return(out)
        } else {
          return(data.frame())
        }
      }, df, c("pos", "neg"), SIMPLIFY = F)
      outlist <- do.call(rbind, outlist)
    } else {
      outlist <- data.frame()
    }
    return(outlist)
  })
  names(modlist) <- compname
  return(modlist)
})
names(countdf_unfilt) <- pathnames

countdf <- countdf_unfilt
modpersub <- lapply(countdf, \(x){
  out <- lapply(x, \(y){
    unique(y$mod)
  })
})
modpersub <- mapply(\(x, y){
  intersect(x, y)
}, modpersub[[1]], modpersub[[2]], SIMPLIFY = F)

countdf <- lapply(countdf, \(x){
  out <- mapply(\(x, y, z){
    out1 <- x |> 
      mutate("comp" = y) |>
      dplyr::filter(mod %in% z)
    return(out1)
  }, x, 
    names(x), modpersub, SIMPLIFY = F)
  return(out)
}) 

countdfjoin <- mapply(\(x, y){
  inner_join(x, y, by = join_by("mod" == "mod"), relationship = "many-to-many") |>
    dplyr::filter(type.x == type.y)
}, countdf[[1]], countdf[[2]], SIMPLIFY = F)

## total # of mods that are:
# - shared between CMC and SZBD
# - match direction (pos/neg)
unlist(lapply(countdfjoin, \(x) length(unique(x$mod))))
  #  conVSCZ
  #      54

# How many celltypes on average are significant among the above mods?
cts_per_mod <- lapply(countdfjoin, \(z){
  temp <- tapply(z$subclass.x, z$mod, list)
  return(lapply(temp, \(x) unique(x)))
})
lapply(cts_per_mod, \(x){
  lapply(x, length) |> unlist() |> summary()
})
# $conVsSCZ
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.000   1.000   1.000   1.796   2.000   7.000 

# Which celltypes are implicated?
mods_per_ct <- lapply(countdfjoin, \(z){
  zpos <- z |> dplyr::filter(type.x=="pos")
  zneg <- z |> dplyr::filter(type.x=="neg")

  temppos <- tapply(zpos$mod, zpos$subclass.x,  list)
  tempneg <- tapply(zneg$mod, zneg$subclass.x,  list)
  return(list("pos" = lapply(temppos, unique),
              "neg" = lapply(tempneg, unique)))
})
lapply(mods_per_ct, \(x){
  lapply(x, \(y) lapply(y, length) |> unlist() |> sort(decreasing = T))
})
# $conVsSCZ
# $conVsSCZ$pos
#      Pvalb      Astro        OPC      L6 CT      L6 IT Chandelier        L6b 
#         11          6          6          5          4          3          3 
#      Micro     Immune L6 IT Car3      Lamp5 Lamp5 Lhx6       Sncg    L2/3 IT 
#          3          2          2          2          2          2          1 
#      L5 IT       Pax6         PC  Sst Chodl 
#          1          1          1          1 

# $conVsSCZ$neg
#     L4 IT     Oligo     L5 IT   L5/6 NP   L2/3 IT      Pax6       Sst Sst Chodl 
#        19         9         3         3         2         2         1         1 
#       Vip 
#         1 

qsave(countdfjoin, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_cmc_szbd_dcopa_overlap.qs"))

paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_CMC"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_SZBDMultiseq"))
pathnames <- c("CMC", "SZBD")

find_overlapping_mods(paths = paths,
                      pathnames = pathnames,
                      common_pool = NULL,
                      save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_cmc_szbd_dcopa_overlap.qs"))
# Fisher p-value:
# [1] 0.934618
 
############
# Run pipeline for overlap (CMC, SZBD)
##############
# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

overlap1 <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_cmc_szbd_dcopa_overlap.qs"))

save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/brainscope_CMC_SZBD_overlap/conVsSCZ")
copa_dir1 <- list(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/brainSCOPE_CMC/"),
                  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/brainSCOPE_SZBDMultiseq/"))
over <- tapply(overlap1[[1]]$mod, overlap1[[1]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("SCZ", "Con"),
               subset = over,
               header_names = c("brainSCOPE_CMC", "brainSCOPE_SZBDMultiseq"))
 



##########################
# SCZ (Brainseq SCZ modules onto all brainSCOPE data)
##########################
 
# Load bulk megaset modules and expr
bulk_expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SCZ_SampleNetworks/1_04-04-29/brainseq_samp_filt_SCZ_1_171_outliers_removed_geneSymbolsAdded.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/")
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/scz/Brainseq_scz")

objs2 <- find_cons_mods(save_path = save_path,
               copa_dir1 = copa_dir1,
               sn_summary_path = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/",
               return_obj = T)



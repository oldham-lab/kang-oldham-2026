# Table structure:
#              Type:    |      Ctrl vs all AD  |  Control vs Early                |  Early vs Late
#                           Up in AD | down in AD:   Up in Early | down in Early:        Up in Late | down in Late:
# Dataset:
# SEAAD2024                     10          2            36               2                0              14
# Morabito [other AD dataset]   3          5
# Rosmap AD [bulk AD mods]      14         19   

library(qs)
library(tidyverse)
library(gt)
library(gtExtras)
library(eulerr)
library(webshot2) # saving gt tables as image
library(showtext)
showtext_auto()

#############################
# Load and count deCOPA (AD)
##############################

paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/rosmapAD_unnormalized"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_rosmap"))
pathnames <- c("bulk_megaset", "rosmapAD", "bulk_MIT", "rosmap_MIT")
 
# List of significant modules
pathList <- lapply(paths, \(d){
  compname <- list.files(d)
  modlist <- lapply(compname, \(x){
    path <- file.path(d, x, "modlist_filtered.qs")
    if(file.exists(path)){
      df <- qread(file.path(d, x, "modlist_filtered.qs"))
      outlist <- lapply(df, \(a) a$mod |> unique())
      names(outlist) <- c("pos", "neg")
    } else {
      outlist <- list()
    }
    return(outlist)
  })
  names(modlist) <- compname
  return(modlist)
})
names(pathList) <- pathnames

# Counts of significant modules
pathCounts <- lapply(paths, \(d){
  compname <- list.files(d)
  modlist <- lapply(compname, \(x){
    path <- file.path(d, x, "modlist_filtered.qs")
    if(file.exists(path)){
      df <- qread(path)
      outvec <- unlist(lapply(df, nrow))
      names(outvec) <- c("pos", "neg")
    } else {
      outvec <- c()
    }
    return(outvec)
  })
  names(modlist) <- compname
  return(modlist)
})
names(pathCounts) <- pathnames

# List of significant modules by subclass
pathList_sub <- lapply(paths, \(d){
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
names(pathList_sub) <- pathnames

###########
# Table of significant mods per subclass
###########
countdf <- pathList_sub[c(1,3)]
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
  }, x[-3], # remove conVLate
    names(x)[-3], modpersub[-3], SIMPLIFY = F)
  return(out)
}) 

countdfjoin <- mapply(\(x, y){
  inner_join(x, y, by = join_by("mod" == "mod")) |>
    dplyr::filter(type.x == type.y)
}, countdf[[1]], countdf[[2]], SIMPLIFY = F)

## total # of mods that are:
# - shared between SEA and MIT
# - match direction (pos/neg)
unlist(lapply(countdfjoin, \(x) length(unique(x$mod))))
  #  conVAll  conVEarly earlyVLate 
  #      158        116          7 

# How many celltypes on average are significant among the above mods?
cts_per_mod <- lapply(countdfjoin, \(z){
  temp <- tapply(z$subclass.x, z$mod, list)
  return(lapply(temp, \(x) unique(x)))
})
lapply(cts_per_mod, \(x){
  lapply(x, length) |> unlist() |> summary()
})
# $conVAll
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.000   2.000   5.000   5.019   7.000  15.000 

# $conVEarly
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.000   1.000   2.000   3.121   4.000  15.000 

# $earlyVLate
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.000   1.000   1.000   1.429   1.500   3.000 

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

# $conVAll        
# $conVAll$pos    
# Endothelial 
#           2                                                                                                                                                            
# $conVAll$neg                                                                                                                                               
#           L4 IT         L5/6 NP           Pvalb             Sst           L5 IT 
#             108              96              93              75              72 
#      L6 IT Car3           L5 ET           Lamp5            Pax6      Lamp5 Lhx6 
#              65              61              60              34              25 
#             Vip           L6 IT            Sncg           L6 CT         L2/3 IT 
#              22              17              14              13              12 
#             L6b             OPC       Astrocyte       Sst Chodl     Endothelial 
#               8               6               4               4               1 
# Oligodendrocyte 
#               1   
# $conVEarly               
# $conVEarly$pos
#           L4 IT           Pvalb             Sst         L5/6 NP           L5 ET 
#              75              41              41              39              30 
#           L5 IT      L6 IT Car3           Lamp5            Pax6           L6 CT 
#              25              24              18              10               9 
#             Vip           L6 IT             OPC            Sncg       Sst Chodl 
#               8               5               4               4               4 
#       Astrocyte         L2/3 IT            VLMC             L6b      Lamp5 Lhx6 
#               3               3               3               2               1 
# Oligodendrocyte 
#               1 
# $conVEarly$neg
#      Chandelier           L6 CT           L6 IT     Endothelial         L2/3 IT 
#               2               2               2               1               1 
#             L6b      Lamp5 Lhx6 Oligodendrocyte            Sncg 
#               1               1               1               1 
# $earlyVLate
# $earlyVLate$pos
#       L6 CT        VLMC   Astrocyte Endothelial         OPC       Pvalb 
#           2           2           1           1           1           1 

# $earlyVLate$neg
#  L6b Pax6 
#    1    1 


# #sea celltypes:
# lapply(countdfjoin, \(x) x$subclass.x) |> unlist() |> unique()
#  [1] "Endothelial" "L4 IT"       "L5 ET"       "L5 IT"       "L5/6 NP"    
#  [6] "Lamp5"       "Pvalb"       "Sst"         "L6 IT"       "L6 IT Car3" 
# [11] "L6 CT"       "Pax6"        "Astrocyte"   "OPC"         "Lamp5 Lhx6" 
# [16] "Sncg"        "Vip"         "L6b"         "Sst Chodl"   "L2/3 IT"    
# [21] "VLMC"        "Chandelier" 
# # mit celltypes:
# lapply(countdfjoin, \(x) x$subclass.y) |> unlist() |> unique()
#  [1] "SMC"              "End"              "Inh VIP"          "Exc L4-5 IT-2"   
#  [5] "Exc L6 IT"        "Exc L3-4 IT"      "Exc L4-5 IT-1"    "Exc L5/6 IT Car3"
#  [9] "Exc L5 ET"        "Inh PVALB"        "Exc L3-5 IT"      "Inh LAMP5"       
# [13] "Exc L6b"          "Exc L2-3 IT"      "OPC"              "Inh PAX6"        
# [17] "Exc L5/6 NP"      "Ast"              "Exc EC"           "Oli"             
# [21] "Exc L5-6 IT"      "Inh SST"          "Exc L6 CT"        "Per" 
common_pool <- list(
  c("Endothelial", "SMC", "VLMC", "End", "Per"),
  c("L4 IT", "Exc L4-5 IT-2", "Exc L3-4 IT","Exc L4-5 IT-1"),
  c("L5 ET",  "Exc L5 ET"),
  c("L5 IT", "Exc L4-5 IT-2", "Exc L4-5 IT-1","Exc L3-5 IT", "Exc L5-6 IT"),
  c("L5/6 NP", "Exc L5/6 NP"),
  c("Lamp5", "Inh LAMP5"),
  c("Pvalb", "Inh PVALB"),
  c("Sst", "Inh SST"),
  c("L6 IT", "Exc L5-6 IT"),
  c("L6 IT Car3", "Exc L5/6 IT Car3"),
  c("L6 CT", "Exc L6 CT"),
  c("Pax6","Inh PAX6"),
  c("Astrocyte", "Ast"),
  c("OPC", "OPC"),
  c("Vip",  "Inh VIP"),
  c("L6b", "Exc L6b"),
  c("L2/3 IT", "Exc L2-3 IT")
)
countdfjoinfilt <- lapply(countdfjoin, \(x){
  these <- lapply(common_pool, \(y){
    which(x$subclass.x %in% y & x$subclass.y %in% y)
  }) |> unlist()
  return(x[these, ])
})

lapply(countdfjoinfilt, \(x) length(unique(x$mod)))
# $conVAll
# [1] 99

# $conVEarly
# [1] 31

# $earlyVLate
# [1] 1

qsave(countdfjoinfilt, file = "seaad2024_mit_dcopa_overlap.qs")

# Which celltypes are implicated?
mods_per_ct_filt <- lapply(countdfjoinfilt, \(z){
  zpos <- z |> dplyr::filter(type.x=="pos")
  zneg <- z |> dplyr::filter(type.x=="neg")

  temppos <- tapply(zpos$mod, zpos$subclass.x,  list)
  tempneg <- tapply(zneg$mod, zneg$subclass.x,  list)
  return(list("pos" = lapply(temppos, unique),
              "neg" = lapply(tempneg, unique)))
})
lapply(mods_per_ct_filt, \(x){
  lapply(x, \(y) lapply(y, length) |> unlist() |> sort(decreasing = T))
})
# $conVAll$pos
# Endothelial    
#           2       
# $conVAll$neg
#      L4 IT      L5 IT L6 IT Car3      Pvalb      L5 ET      Lamp5       Pax6 
#         48         41         37         23         19         19         14 
#        Vip    L5/6 NP    L2/3 IT 
#          6          5          2 
# $conVEarly
# $conVEarly$pos
#      L5 ET L6 IT Car3      L4 IT      L5 IT        OPC       Pax6 
#         15         12          7          5          2          2 
# $conVEarly$neg
# Endothelial 
#           1 
# $earlyVLate
# $earlyVLate$pos
# L6 CT 
#     1 
# $earlyVLate$neg
# NULL

####
# AD table
####

ad_tab <- data.frame(
  "Dataset" = c("SEAAD2024 (Bulk megaset)", "SEAAD2024 (ROSMAP)", "MIT (Bulk megaset)", "MIT (ROSMAP)", "Morabito 2021 (Bulk megaset)"),
  "Up in AD vs con" = c(9, 9, 2, 0, 3), 
  "Down in AD vs con" = c(42, 44, 89, 83, 4),
  "Up in early AD vs con" = c(4, 0, 3, 1, NA),
  "Down in early AD vs con" = c(52, 44, 102, 89, NA),
  "Up in late AD vs con" = c(7, 12, 2, 1, NA),
  "Down in late AD vs con" = c(26, 31, 30, 34, NA),
  "Up in late AD vs early" = c(4, 8, 5, 2, NA), 
  "Down in late AD vs early" = c(4, 7, 13, 12, NA)
)
 
tab1 <- gt(ad_tab, rowname_col="genes", groupname_col="celltype") |>
  tab_spanner(
    label = "Con vs all AD",
    columns = c("Up.in.AD.vs.con", "Down.in.AD.vs.con")
  ) |>
  tab_spanner(
    label = "Con vs early AD",
    columns = c("Up.in.early.AD.vs.con", "Down.in.early.AD.vs.con")
  ) |>
  tab_spanner(
    label = "Con vs late AD",
    columns = c("Up.in.late.AD.vs.con", "Down.in.late.AD.vs.con")
  ) |>
  tab_spanner(
    label = "Early vs late",
    columns = c("Up.in.late.AD.vs.early", "Down.in.late.AD.vs.early")
  ) |>
  cols_label(
    Up.in.AD.vs.con = html("Up in AD vs con"),
    Down.in.AD.vs.con = html("Down in AD vs con"),
    Up.in.early.AD.vs.con = html("Up in early AD vs con"),
    Down.in.early.AD.vs.con = html("Down in early AD vs con"),
    Up.in.late.AD.vs.con = html("Up in late AD vs con"),
    Down.in.late.AD.vs.con = html("Down in late AD vs con"),
    Up.in.late.AD.vs.early = html("Up in early AD"),
    Down.in.late.AD.vs.early = html("Down in early AD")
  )

gtsave(tab1, filename = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/summary_table.html"))

#######
# Count deCOPA mods before filtering
#######
# SEAAD2024 (con vs all)
copa_dir_sea <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
sea_euc_dist <- list("all"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_all.qs")),
                     "pos"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                     "neg"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_negative.qs")))
a1 <- sea_euc_dist$all[[3]]$all
# > length(a1)
# [1] 217
s1 <- sea_euc_dist$all[[3]][which(names(sea_euc_dist$all[[3]]) != "all")] |> unlist() |> unique()
# > sum(a1 %in% s1)
# [1] 217
# > sum(!s1 %in% a1)
# [1] 161
# 217 modules show significant differences in module projections over all subclasses
# 161 modules show significant differences in at least one individual subclass, but not across all subclasses
sr1 <- sea_euc_dist$all[[3]][which(names(sea_euc_dist$all[[3]]) != "all")] |>
  lapply(\(x) x[!x %in% a1])
countvec <- c()
for(i in 1:length(sr1)){
  q1 <- unlist(sr1[-i])
  countvec[i] <- sum(!sr1[[i]] %in% q1)
}
# > sum(countvec)
# [1] 101
# Of those 161 modules, 101 only show significant differences in a single subclass

# MIT (con vs all)
copa_dir_mit <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/")
mit_euc_dist <- list("all"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_all.qs")),
                     "pos"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                     "neg"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_negative.qs")))
am1 <- mit_euc_dist$all[[3]]$all
# length(am1)
# [1] 176
sm1 <- mit_euc_dist$all[[3]][which(names(mit_euc_dist$all[[3]]) != "all")] |> unlist() |> unique()
# sum(am1 %in% sm1)
# [1] 176
# sum(!sm1 %in% am1)
# [1] 215
# 176 modules show significant differences in module projections over all subclasses
# 215 modules show significant differences in at least one individual subclass, but not across all subclasses
srm1 <- mit_euc_dist$all[[3]][which(names(mit_euc_dist$all[[3]]) != "all")] |>
  lapply(\(x) x[!x %in% am1])
countvecm <- c()
for(i in 1:length(srm1)){
  q1 <- unlist(srm1[-i])
  countvecm[i] <- sum(!srm1[[i]] %in% q1)
}
# > sum(countvecm)
# [1] 123
# Of those 215 modules, 123 show significant differences only in a single subclass

# Comparing SEA vs MIT
ai1 <- intersect(a1, am1)
# length(ai1)
# > 139
si1 <- intersect(s1, sm1)
complist <- list()
complist[[1]] <- si1

# length(si1)
# > 288
# sum(ai1 %in% si1)
# > 139 
# 288 modules show a significant difference between control and AD.
# (Of those, 139 are significantly different when considering all celltypes and for at least one individual subclass)
# (288-139 = 149 of those are not significant across all celltypes but are for individual subclasses)

decopaov <- intersect(pathList$bulk_megaset$conVAll |> unlist(), pathList$bulk_MIT$conVAll |> unlist())
# > sum(decopaov %in% si1)
# [1] 31
# > length(decopaov)
# [1] 31
# As expected all of the filtered modules are part of the unfiltered significant modules 

## Con vs early
copa_dir_sea <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/")
sea_euc_dist <- list("all"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_all.qs")),
                     "pos"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                     "neg"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_negative.qs")))
a1 <- sea_euc_dist$all[[3]]$all # 228 mods
s1 <- sea_euc_dist$all[[3]][which(names(sea_euc_dist$all[[3]]) != "all")] |> unlist() |> unique() # 382 mods
# 224 modules show significant differences in module projections over all subclasses and in individual subclasses
# 154 modules show significant differences in at least one individual subclass, but not across all subclasses
# 4 modules show significant differences in module projections over all subclasses but not in individual subclasses
sr1 <- sea_euc_dist$all[[3]][which(names(sea_euc_dist$all[[3]]) != "all")] |>
  lapply(\(x) x[!x %in% a1])
countvec <- c()
for(i in 1:length(sr1)){
  q1 <- unlist(sr1[-i])
  countvec[i] <- sum(!sr1[[i]] %in% q1)
}
# Of those 158 modules, 112 only show significant differences in a single subclass
# MIT 
copa_dir_mit <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly/")
mit_euc_dist <- list("all"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_all.qs")),
                     "pos"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                     "neg"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_negative.qs")))
am1 <- mit_euc_dist$all[[3]]$all # 199 mods
sm1 <- mit_euc_dist$all[[3]][which(names(mit_euc_dist$all[[3]]) != "all")] |> unlist() |> unique() # 415 mods
ai1 <- intersect(a1, am1) # 148 mods
si1 <- intersect(s1, sm1) # 286 mods
complist[[2]] <- si1
decopaov <- intersect(pathList$bulk_megaset$conVEarly |> unlist(), pathList$bulk_MIT$conVEarly |> unlist())
sum(decopaov %in% si1) # 21 mods

## Early vs late
copa_dir_sea <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/")
sea_euc_dist <- list("all"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_all.qs")),
                     "pos"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                     "neg"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_negative.qs")))
a1 <- sea_euc_dist$all[[3]]$all # 80 mods
s1 <- sea_euc_dist$all[[3]][which(names(sea_euc_dist$all[[3]]) != "all")] |> unlist() |> unique() # 209 mods
# 79 modules show significant differences in module projections over all subclasses and in at least one individual subclass
# 130 modules show significant differences in at least one individual subclass, but not across all subclasses
# 1 module shows significant differences in module projections over all subclasses but not in at least one individual subclass
sr1 <- sea_euc_dist$all[[3]][which(names(sea_euc_dist$all[[3]]) != "all")] |>
  lapply(\(x) x[!x %in% a1])
countvec <- c()
for(i in 1:length(sr1)){
  q1 <- unlist(sr1[-i])
  countvec[i] <- sum(!sr1[[i]] %in% q1)
}
# Of those 158 modules, 112 only show significant differences in a single subclass
# MIT 
copa_dir_mit <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate/")
mit_euc_dist <- list("all"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_all.qs")),
                     "pos"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                     "neg"=qread(paste0(copa_dir_mit, "/euclidean_distances/euclidean_sigmods_negative.qs")))
am1 <- mit_euc_dist$all[[3]]$all # 162 mods
sm1 <- mit_euc_dist$all[[3]][which(names(mit_euc_dist$all[[3]]) != "all")] |> unlist() |> unique() # 332 mods
ai1 <- intersect(a1, am1) # 42 mods
si1 <- intersect(s1, sm1) # 135 mods
complist[[3]] <- si1
decopaov <- intersect(pathList$bulk_megaset$earlyVLate |> unlist(), pathList$bulk_MIT$earlyVLate |> unlist())
sum(decopaov %in% si1) # 1 mod

## Overlap between con/ad, con/earlyad, early/latead
length(intersect(complist[[1]], complist[[2]]))
# 241
length(intersect(complist[[1]], complist[[3]]))
# 104
length(intersect(complist[[2]], complist[[3]]))
# 104

############
# Find # of modules where significant mod is same as highest mod
##########
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)
# SEA
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
euc_dist <- list("all"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(paste0(copa_dir_sea, "/euclidean_distances/euclidean_sigmods_negative.qs")))
sea_means <- qread(paste0(copa_dir1, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
sea_means[[2]] <- sea_means[[2]][match(rownames(sea_means[[1]]), rownames(sea_means[[2]])), ]
sea_means <- lapply(sea_means, as.data.frame)
sub_diff <- sea_means[[1]] - sea_means[[2]] |> as.data.frame()
seapval <- fread(data.table = F, file = paste0(copa_dir1, "/euclidean_distances/p_values_all_Subclass.csv"))
seadist <- fread(data.table = F, file = paste0(copa_dir1, "/euclidean_distances/euclidean_distances_all_Subclass.csv"))
rownames(seadist) <- these_mods
rownames(seapval) <- these_mods
# all mods where max expression subclass matches decopa significant subclass
copamaxallsea <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(modind){ # Add info about subclass with highest expression in AD
              mean1 <- sea_means[[1]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
              return(names(mean1)[which.max(mean1)])
            }) |> unlist() 
  if(length(outvec) > 0){
    return(data.frame("subclass" = subclass_name, "maxct" = outvec))
  } else {
    return(data.frame("subclass" = subclass_name, "maxct" = NA))
  }
}, euc_dist$all[[3]], names(euc_dist$all[[3]]), SIMPLIFY = F) 
copamaxmatchcountsea <- lapply(copamaxallsea, \(x) rownames(x)[x$subclass == x$maxct])
sum(lapply(copamaxmatchcountsea, length) |>unlist())
# [1] 202
# all mods with coherent mod genes (regardless of direction)
copacohallsea <- mapply(\(subclass, subclass_name){
    outvec <- lapply(subclass, \(mod){
        temp <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], colnames(sub_diff) == subclass_name]
        return(sum(temp > 0) / length(temp))
    }) |> unlist() 
    if(length(outvec) > 0){
      return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
    } else {
      return(data.frame("subclass" = subclass_name, "pcnt" = NA))
    }
}, euc_dist$all[[3]], names(euc_dist$all[[3]]), SIMPLIFY = F) |>
  lapply(\(x){
    x[(x$pcnt == 1) | (x$pcnt == 0), ] |>
    rownames_to_column(var = "mod")
  }) |> do.call(what = "rbind")  |>
  dplyr::filter(!is.na(pcnt))

#sum(copamaxmatchcount |> unlist() %in% copacohall$mod)
# 163

# MIT
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/")
euc_dist <- list("all"=qread(paste0(copa_dir1, "/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(paste0(copa_dir1, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(paste0(copa_dir1, "/euclidean_distances/euclidean_sigmods_negative.qs")))
sea_means <- qread(paste0(copa_dir1, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
sea_means[[2]] <- sea_means[[2]][match(rownames(sea_means[[1]]), rownames(sea_means[[2]])), ]
sea_means <- lapply(sea_means, as.data.frame)
sub_diff <- sea_means[[1]] - sea_means[[2]] |> as.data.frame()
seapval <- fread(data.table = F, file = paste0(copa_dir1, "/euclidean_distances/p_values_all_Subclass.csv"))
seadist <- fread(data.table = F, file = paste0(copa_dir1, "/euclidean_distances/euclidean_distances_all_Subclass.csv"))
rownames(seadist) <- these_mods
rownames(seapval) <- these_mods
# all mods where max expression subclass matches decopa significant subclass
copamaxallmit <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(modind){ # Add info about subclass with highest expression in AD
              mean1 <- sea_means[[1]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
              return(names(mean1)[which.max(mean1)])
            }) |> unlist() 
  if(length(outvec) > 0){
    return(data.frame("subclass" = subclass_name, "maxct" = outvec))
  } else {
    return(data.frame("subclass" = subclass_name, "maxct" = NA))
  }
}, euc_dist$all[[3]], names(euc_dist$all[[3]]), SIMPLIFY = F) 
cmmcountmit <- lapply(copamaxallmit, \(x) rownames(x)[x$subclass == x$maxct])
sum(lapply(cmmcountmit, length) |>unlist())
# [1] 209
# all mods with coherent mod genes (regardless of direction)
copacohallmit <- mapply(\(subclass, subclass_name){
    outvec <- lapply(subclass, \(mod){
        temp <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], colnames(sub_diff) == subclass_name]
        return(sum(temp > 0) / length(temp))
    }) |> unlist() 
    if(length(outvec) > 0){
      return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
    } else {
      return(data.frame("subclass" = subclass_name, "pcnt" = NA))
    }
}, euc_dist$all[[3]], names(euc_dist$all[[3]]), SIMPLIFY = F) |>
  lapply(\(x){
    x[(x$pcnt == 1) | (x$pcnt == 0), ] |>
    rownames_to_column(var = "mod")
  }) |> do.call(what = "rbind") |>
  dplyr::filter(!is.na(copacohallmit))

maxoverlap <- intersect(copamaxmatchcountsea |> unlist(), cmmcountmit |> unlist())
cohoverlap <- intersect(copacohallsea$mod, copacohallmit$mod)
alleulermods <- unique(c(maxoverlap, cohoverlap, complist[[1]]))

eulermat <- data.frame("Significantly different" = alleulermods %in% complist[[1]],
                       "Highest expressed subclass is sig" = alleulermods %in% maxoverlap,
                       "Consistent module gene change" = alleulermods %in% cohoverlap)
fit <- eulerr::euler(eulermat)
pdf(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/eulerr/convall_overview.pdf")))
p <- plot(fit,
          quantities = list(type = c("counts", "percent"), cex = 2),
          labels = list(font = 1, cex = 4))
print(p)
dev.off()

pdf(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/eulerr/convall_overview.pdf_noLabs.pdf")))
p2 <- plot(fit,
           fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
           quantities = list(type = c("counts", "percent"), cex = 2),
           labels = F,
           edges = F,
           family = "sans")
print(p2)
dev.off()



####
# SCZ
####

scz_tab <- data.frame(
  "Dataset" = c("brainSCOPE", "brainSCOPE (Brainseq SCZ modules)"),
  "Up in SCZ" = c(3, 7), 
  "Down in SCZ" = c(39, 44)
)

tab2 <- gt(scz_tab, rowname_col="genes", groupname_col="celltype") |>
  cols_label(
    Up.in.SCZ = html("Up in SCZ"),
    Down.in.SCZ = html("Down in SCZ"),
  )

gtsave(tab2, filename = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/summary_table_scz.html"))


######
# Overlap between SEAAD2024 and MIT
#################

# Compare bulk megaset modules
overlap_list <- list()
overlap_pctg <- list()
for(i in 1:4){
  outlist <- list()
  outlist[[1]] <- intersect(pathList$bulk_megaset[[i]][[1]], pathList$bulk_MIT[[i]][[1]])
  outlist[[2]] <- intersect(pathList$bulk_megaset[[i]][[2]], pathList$bulk_MIT[[i]][[2]])
  names(outlist) <- c("pos", "neg")
  overlap_list[[i]] <- outlist

  t1 <- length(unique(c(pathList$bulk_megaset[[i]][[1]], pathList$bulk_MIT[[i]][[1]])))
  t2 <- length(unique(c(pathList$bulk_megaset[[i]][[2]], pathList$bulk_MIT[[i]][[2]])))

  overlap_pctg[[i]] <- data.frame("comparison" = names(pathList$bulk_megaset)[i],
                                  "direction" = c("pos", "neg"), 
                                  "count" = c(length(outlist[[1]]), length(outlist[[2]])),
                                  "pctg" = c(length(outlist[[1]])/t1, length(outlist[[2]])/t2))

}
names(overlap_pctg) <- names(pathList$bulk_megaset)
overlap_pctg <- do.call(rbind, overlap_pctg)

p <- ggplot(overlap_pctg, aes(x = comparison, y = count)) +
  theme_classic() + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Overlapping\nmodules") 
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/overlap_bulkMegaset_ADdatasets.svg"), bg = "white")
 

# Compare bulk megaset modules using grouped barchart
comparison_names <- c("Con vs all AD", "Con vs early AD", "Con vs Late AD", "Early vs Late AD")
#overlap_list <- list()
overlap_pctg <- list()
for(i in 1:4){
  outlist <- list()
  outlist[[1]] <- intersect(pathList$bulk_megaset[[i]][[1]], pathList$bulk_MIT[[i]][[1]])
  outlist[[2]] <- intersect(pathList$bulk_megaset[[i]][[2]], pathList$bulk_MIT[[i]][[2]])
  names(outlist) <- c("pos", "neg")
  #overlap_list[[i]] <- outlist

  t1 <- length(unique(c(pathList$bulk_megaset[[i]][[1]], pathList$bulk_MIT[[i]][[1]])))
  t2 <- length(unique(c(pathList$bulk_megaset[[i]][[2]], pathList$bulk_MIT[[i]][[2]])))

  overlap_pctg[[i]] <- data.frame("comparison" = comparison_names[i],
                                 # "direction" = c("pos", "neg"), 
                                  "SEAAD2024" = length(pathList$bulk_megaset[[i]][[1]]) + length(pathList$bulk_megaset[[i]][[2]]),
                                  "MIT" = length(pathList$bulk_MIT[[i]][[1]]) + length(pathList$bulk_MIT[[i]][[2]]),
                                  "shared" = length(outlist[[1]]) + length(outlist[[2]]))

}
names(overlap_pctg) <- names(pathList$bulk_megaset)
overlap_pctg <- do.call(rbind, overlap_pctg)
overlap_pctg <- pivot_longer(overlap_pctg, !comparison, names_to = "which", values_to = "count")
overlap_pctg$which <- factor(overlap_pctg$which, levels = unique(overlap_pctg$which))

p <- ggplot(overlap_pctg, aes(x = comparison, y = count, fill = which)) +
  theme_classic() + 
  geom_bar(stat = "identity", position = "stack", alpha = 0.6) +
  labs(x = "", y = "# of significant modules", title = "Bulk megaset modules") +
  theme(legend.title = element_blank(), 
        text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 20, angle = 30, hjust = 1, vjust = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/overlap_bulkMegaset_ADdatasets_stacked.svg"), bg = "white")
 
# Create venn diagrams for each comparison
for(i in 1:4){
  outlist <- list()
  outlist[[1]] <- intersect(pathList$bulk_megaset[[i]][[1]], pathList$bulk_MIT[[i]][[1]])
  outlist[[2]] <- intersect(pathList$bulk_megaset[[i]][[2]], pathList$bulk_MIT[[i]][[2]])
  names(outlist) <- c("pos", "neg")
  #overlap_list[[i]] <- outlist

  t1 <- unique(c(pathList$bulk_megaset[[i]][[1]], pathList$bulk_MIT[[i]][[1]]))
  t2 <- unique(c(pathList$bulk_megaset[[i]][[2]], pathList$bulk_MIT[[i]][[2]]))

  # Create overlap mat for eulerr (combining pos and neg)
  eulerdf <- data.frame("SEAAD2024" = unique(c(t1, t2)) %in% unlist(pathList$bulk_megaset[[i]]),
                        "MIT_Multiome" = unique(c(t1, t2)) %in% unlist(pathList$bulk_MIT[[i]]))
  
  fit <- euler(eulerdf)
  pdf(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/eulerr/"), paste0(names(pathList[[1]])[i], ".pdf")))
  p <- plot(fit,
     quantities = list(type = c("counts", "percent"), cex = 2),
     labels = list(font = 1, cex = 4))
  print(p)
  dev.off()

  pdf(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/eulerr/"), paste0(names(pathList[[1]])[i], "_noLabs.pdf")))
  p2 <- plot(fit,
       fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
       quantities = list(type = c("counts", "percent"), cex = 2),
       labels = F,
       edges = F,
       family = "sans")
  print(p2)
  dev.off()
}

# Create legend
# #0c7175
colframe <- data.frame("cols" = c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"),
                      "vals" = c(1,2,3),
                      "names" = c("SEAAD2024", "MIT_Multiome", "Both"))                      
colvec <- colframe$cols
names(colvec) <- colframe$names
colframe$names <- factor(colframe$names, levels = colframe$names)
p <- ggplot(colframe, aes(x = cols, y = vals, fill = names)) + 
  geom_bar(stat='identity', alpha = 0.6) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = colvec)
pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/eulerr/legend.pdf"))
legend_grob <- cowplot::get_legend(p)
grid::grid.newpage()
grid::grid.draw(legend_grob)
dev.off()

# Compare ROSMAP modules using grouped barchart
comparison_names <- c("Con vs all AD", "Con vs early AD", "Con vs Late AD", "Early vs Late AD")
#overlap_list <- list()
overlap_pctg <- list()
for(i in 1:4){
  outlist <- list()
  outlist[[1]] <- intersect(pathList$rosmapAD[[i]][[1]], pathList$rosmap_MIT[[i]][[1]])
  outlist[[2]] <- intersect(pathList$rosmapAD[[i]][[2]], pathList$rosmap_MIT[[i]][[2]])
  names(outlist) <- c("pos", "neg")
  #overlap_list[[i]] <- outlist

  t1 <- length(unique(c(pathList$rosmapAD[[i]][[1]], pathList$rosmap_MIT[[i]][[1]])))
  t2 <- length(unique(c(pathList$rosmapAD[[i]][[2]], pathList$rosmap_MIT[[i]][[2]])))

  overlap_pctg[[i]] <- data.frame("comparison" = comparison_names[i],
                                 # "direction" = c("pos", "neg"), 
                                  "SEAAD2024" = length(pathList$rosmapAD[[i]][[1]]) + length(pathList$rosmapAD[[i]][[2]]),
                                  "MIT" = length(pathList$rosmap_MIT[[i]][[1]]) + length(pathList$rosmap_MIT[[i]][[2]]),
                                  "shared" = length(outlist[[1]]) + length(outlist[[2]]))

}
names(overlap_pctg) <- names(pathList$rosmapAD)
overlap_pctg <- do.call(rbind, overlap_pctg)
overlap_pctg <- pivot_longer(overlap_pctg, !comparison, names_to = "which", values_to = "count")
overlap_pctg$which <- factor(overlap_pctg$which, levels = unique(overlap_pctg$which))

p <- ggplot(overlap_pctg, aes(x = comparison, y = count, fill = which)) +
  theme_classic() + 
  geom_bar(stat = "identity", position = "stack", alpha = 0.6) +
  labs(x = "", y = "# of significant modules", title = "ROSMAP AD modules") +
  theme(legend.title = element_blank(), 
        text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 20, angle = 30, hjust = 1, vjust = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/overlap_rosmap_ADdatasets_stacked.svg"), bg = "white")


# Figure outline:
# a) schematic
# b) table of overlaps
# c) flowchart of celltypes (con early late)
# d) examples
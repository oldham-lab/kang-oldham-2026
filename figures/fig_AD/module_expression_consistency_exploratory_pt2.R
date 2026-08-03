# Of the modules that are expressed differently in control vs disease, what module genes are driving this?
# Code for COPA results specifically

library(DESeq2)
library(qs)
library(data.table)
library(AnnotationHub)
library(tidyverse)
library(cowplot)
options(bitmapType = 'cairo')

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_expression_consistency_fxns.R"))

##########################
# Load objects of interest
##########################

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Plot
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized/conVAll")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T)
 
### What is significance of difference between con and ad over all genes per celltype?
plot_allgenes_diff(copa_dir1 = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/"),
                   ytitle = "Mean expression difference\n(AD minus control)",
                   plottitle = "AD vs con expression over all genes in SEAAD2024",
                   save_dir = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/ADvsCon_over_allGenes.png"))

plot_allgenes_diff(copa_dir1 = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly"),
                   ytitle = "Mean expression difference\n(Con minus early AD)",
                   plottitle = "Con vs early AD expression over all genes in SEAAD2024",
                   save_dir = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/conVsEarly_over_allGenes.png"))

plot_allgenes_diff(copa_dir1 = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate"),
                   ytitle = "Mean expression difference\n(Late minus early AD)",
                   plottitle = "Early vs late AD expression over all genes in SEAAD2024",
                   save_dir = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/earlyVsLate_over_allGenes.png"))

plot_allgenes_diff(copa_dir1 = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/",
                   ytitle = "Mean expression difference\n(SCZ minus con)",
                   plottitle = "Con vs SCZ expression over all genes in brainSCOPE",
                   save_dir = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/conVsSCZ_over_allGenes.png"))

plot_allgenes_diff(copa_dir1 = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/"),
                   ytitle = "Mean expression difference\n(AD minus con)",
                   plottitle = "AD vs Con expression over all genes in MIT AD Multiome/Multiregion (MTC)",
                   save_dir = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/adVSCon_MITMultiomePFC_over_allGenes.png"))
 

# Con vs late
copa_dir1 = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly")
copa_dir2 = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate")
ytitle = "Mean expression difference\n(Con minus late AD)"
plottitle = "Con vs late AD expression over all genes in SEAAD2024"
save_dir = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/conVsLateAD_over_allGenes.png")
sea_means <- qread(paste0(copa_dir1, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
sea_means[[2]] <- sea_means[[2]][match(rownames(sea_means[[1]]), rownames(sea_means[[2]])), ]
sea_means <- lapply(sea_means, as.data.frame)
sea_means2 <- qread(paste0(copa_dir2, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
sea_means2[[2]] <- sea_means2[[2]][match(rownames(sea_means2[[1]]), rownames(sea_means2[[2]])), ]
sea_means2 <- lapply(sea_means2, as.data.frame)
sea_means[[1]] <- sea_means[[1]][, colnames(sea_means[[1]]) %in% colnames(sea_means2[[1]])]

sub_diff <- sea_means[[1]] - sea_means2[[1]] |> as.data.frame()
sub_diff_means <- apply(sub_diff, 2, mean)

meantestsbonf <- mapply(\(ad, con){
  p <- wilcox.test(ad, con)$p.value
  return(p * 22)
}, sea_means[[1]] |> as.list(), sea_means2[[1]] |> as.list(), SIMPLIFY = F) |> unlist()

plotdf <- data.frame("ct" = names(sub_diff_means), "meandiff" = sub_diff_means, "pval" = -log10(meantestsbonf))
plotdf$star <- ""
plotdf$star[plotdf$pval >= -log10(.05)] <- "*"
plotdf$star[plotdf$pval >= -log10(.01)] <- "**"
plotdf$star[plotdf$pval >= -log10(.001)] <- "***"
plotdf$ypos = ifelse(plotdf$meandiff > 0, plotdf$meandiff, 0)

p <- ggplot(plotdf, aes(x = ct, y = meandiff, fill = pval)) + 
  theme_classic() +
  geom_bar(stat = "identity") +
  geom_text(aes(label = star, y = ypos), vjust = -0.5) +
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 16),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.title = element_text(hjust = 0.5, size = 12),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = ytitle,
       fill = bquote(-log[10]~"p-val"),
       title = plottitle)
ggsave(p, file = save_dir, width = 11, height = 6.5)




################################
# Coherence in SEAAD2024 modules
################################

# How many modules per subclass are completely coherent (in SEAAD data)?
# Positive modules (higher in AD):
copacohpos <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], colnames(sub_diff) == subclass_name]
    return(sum(temp > 0) / length(temp))
  }) |> unlist() 
  return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
}, euc_dist$pos[[3]], names(euc_dist$pos[[3]]), SIMPLIFY = F) #|> do.call(what = "rbind")# |>
#   ggplot(aes(x = subclass, y = pcnt)) + 
#     theme_bw() + geom_violin() + coord_flip() + theme(text = element_text(size = 12))
#lapply(copacohpos, \(x) sum(x$pcnt == 1)) |> unlist()
#       Astrocyte      Chandelier     Endothelial         L2/3 IT           L4 IT 
#               4               1               8               1               0 
#           L5 ET           L5 IT         L5/6 NP           L6 CT           L6 IT 
#               1               0               0               0               1 
#      L6 IT Car3             L6b           Lamp5      Lamp5 Lhx6   Microglia-PVM 
#               0               1               0               0               0 
# Oligodendrocyte             OPC            Pax6           Pvalb            Sncg 
#               2               0               0               0               0 
#             Sst       Sst Chodl             Vip            VLMC 
#               1               0               1               0 
# Anywhere from 0-4 modules are completely coherent and higher in con for each subclass.
# How many overlap with mods significant in ROSMAP?
#test1 <- lapply(copacohpos, \(x) rownames(x)[x$pcnt==1]) |> unlist() |> as.numeric()
#sum(test1 %in% which(consis_strict_con))
# 0
# How many of these modules are specific to a subclass
copacohpos_complete <- lapply(copacohpos, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
#tapply(copacohpos_complete[ ,2], copacohpos_complete[ ,1], list)
# $`1055`                                                                                                                                                         
# [1] "Endothelial"                                                                                                                               
# $`13`                                                                                                                                                           
# [1] "Endothelial"                                                                                                                                                                                                                                                                                                             
# $`149`                                                                                                                                                          
# [1] "Endothelial"     "Oligodendrocyte"                                                                                                                                                                                        
# $`181`                                                                          
# [1] "L2/3 IT"                                                                                                                                                                                                                                  
# $`182`                                                                                                                                                          
# [1] "Endothelial"                                                                                                       
# $`348`                                                                          
# [1] "Astrocyte"                                                                                                        
# $`474`                               
# [1] "Astrocyte"
# $`514`
# [1] "Chandelier" "L5 ET"      "L6 IT"      "L6b"        "Sst"       
# [6] "Vip"       
# $`540`
# [1] "Oligodendrocyte"
# $`693`
# [1] "Astrocyte"
# $`716`
# [1] "Endothelial"
# $`882`
# [1] "Endothelial"
# $`947`
# [1] "Endothelial"
# $`948`
# [1] "Endothelial"
# $`970`
# [1] "Astrocyte"
# The majority of modules are only completely coherent in one subclass.
# 4 modules are completely coherent in two subclasses
# 1 module (514) is completely coherent in 6 subclasses - this is a mitochondrial module

#  modules (higher in con):
copacohneg <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], colnames(sub_diff) == subclass_name]
    return(sum(temp < 0) / length(temp))
  }) |> unlist() 
  return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
}, euc_dist$neg[[3]], names(euc_dist$neg[[3]]), SIMPLIFY = F) 
#lapply(copacohneg, \(x) sum(x$pcnt == 1)) |> unlist()
#       Astrocyte      Chandelier     Endothelial         L2/3 IT           L4 IT 
#               5               0               3               5              13 
#           L5 ET           L5 IT         L5/6 NP           L6 CT           L6 IT 
#               0              10              17               2               9 
#      L6 IT Car3             L6b           Lamp5      Lamp5 Lhx6   Microglia-PVM 
#               3               2              24               2               0 
# Oligodendrocyte             OPC            Pax6           Pvalb            Sncg 
#               1               0               2               6               1 
#             Sst       Sst Chodl             Vip            VLMC 
#              20               4               5               0 
# More coherent modules that are higher in control than there are in AD.
# - Around 4-6 on average with max of 23 (Lamp5) are completely coherent and higher in AD.
# How many overlap with mods significant in ROSMAP?
#test1 <- lapply(copacohneg, \(x) rownames(x)[x$pcnt==1]) |> unlist() |> as.numeric()
#sum(test1 %in% which(consis_strict_AD))
# 0 
# How many of these modules are specific to a single subclass?
copacohneg_complete <- lapply(copacohneg, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
#copacohneg_count <- tapply(copacohneg_complete[ ,2], copacohneg_complete[ ,1], list) |>
#  lapply(length) |> unlist()
#sum(copacohneg_count > 1)
# 26
#sum(copacohneg_count == 1)
# 31
# 31 of the 57 completely coherent modules higher in AD are specific to a single subclass
#summary(copacohneg_count)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.000   1.000   1.000   2.351   2.000  16.000 
# mod 149 is an outlier that is coherent in 16 subclasses (all neuronal)
# - this module is a ribosomal module that is higher in AD in all 
#test <- tapply(copacohneg_complete[ ,2], copacohneg_complete[ ,1], list) 
#names(euc_dist$neg[[3]])[!names(euc_dist$neg[[3]]) %in% test$`149`]
# "Astrocyte"       "Chandelier"      "Endothelial"     "L5 ET"          
# [5] "Microglia-PVM"   "Oligodendrocyte" "OPC"             "VLMC"  
# mod 149 is not present in glial celltypes and L5ET
# Other mods present in many subclasses:
# mod 182 (9 subclasses, neuronal): electron transport related
# mod 716 (9 subclasses, neuronal + astro): ribosomal related
# mod 295 (7 subclasses, neuronal + astro): ribosomal related
# mod 25 (7 subclasses, neuronal + astro): ETC related
# mod 416 (6 subclasses, neuronal): not sure
# mod 85 (6 subclasses, neuronal + astro): ribosomal
# mod 131 (5 subclasses, excitatory + VIP): no idea
# mod 441 (5 subclasses, neuronal): tubulins
# mod 451 (5 subclasses, excitatory): no idea (protein/AA modification, some mito-related?)

# For majority of subclasses, coherence in AD modules ranges from 50% to 100%, median roughly 80-90% ish
#summary(mod_bc_lengths[euc_pos_fdr])
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  4.00   12.00   18.00   20.13   26.00   53.00 
# All the AD modules fairly large modules as well

#####################
# Outline of figures:
#####################

### Plot the results of the seattle analysis

# Mods higher in AD:
thesead <- unique(copacohpos_complete[,1]) |> as.numeric()
# Which of the coherent mods align with copa mods unique for 1 subclass
eucposcommon <- unlist(euc_dist$pos[[3]])
eucposcommon <- unique(eucposcommon[duplicated(eucposcommon)])
eucposu <- lapply(euc_dist$pos[[3]], \(x){
  return(x[!x %in% eucposcommon])
})
#sum(thesead %in% unlist(eucposu))
# 12 of the 15 coherent modules are significantly unique for a single celltype
thesead2 <- thesead[thesead %in% unlist(eucposu)]
# Calc percentage diff
sub_diff_per <- ((sea_means[[1]] - sea_means[[2]]) / sea_means[[1]]) 

copacohpos_complete$copa_unique <- F
for(i in 1:nrow(copacohpos_complete)){
  temp <- names(euc_dist$pos[[3]])[unlist(lapply(euc_dist$pos[[3]], \(x) copacohpos_complete$mod[i] %in% x))]
  if(length(temp) == 1){
    copacohpos_complete$copa_unique[i] <- temp
  } else if(length(temp) == 2){
    copacohpos_complete$copa_unique[i] <- paste(temp, collapse = " ")
  } else {
    copacohpos_complete$copa_unique[i] <- "multiple"
  }
}
copacohpos_complete <- copacohpos_complete |>
  mutate(max_AD = lapply(mod, \(modind){ # Add info about subclass with highest expression in AD
    mean1 <- sea_means[[1]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
      return(names(mean1)[which.max(mean1)])
    }) |> unlist()) |>
  mutate(max_con = lapply(mod, \(modind){ # Add info about subclass with highest expression in con
    mean1 <- sea_means[[2]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
      return(names(mean1)[which.max(mean1)])
    }) |> unlist()) |>
  mutate(dist_all = lapply(mod, \(modind){ # Add info about total euclidean distance (over all subclasses) 
    seadist$all[which(these_mods %in% modind)]
    }) |> unlist()) |>
  mutate(dist_per = lapply(mod, \(modind){ # Add info about subclass with highest percentage diff
    tempdf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[modind]], ] |> colMeans() |> unlist()
    return(names(tempdf)[which.max(tempdf)])
    }) |> unlist()) |>
  dplyr::filter((copa_unique == max_AD | copa_unique == max_con) & subclass == copa_unique) |>
  arrange(desc(dist_all))
#rownames(copacohpos_complete) <- paste0(rownames(copacohpos_complete), "_", copacohpos_complete$mod)

# Repeat the above steps for negative (higher in con) modules
copacohneg_complete$copa_unique <- F
for(i in 1:nrow(copacohneg_complete)){
  temp <- names(euc_dist$pos[[3]])[unlist(lapply(euc_dist$pos[[3]], \(x) copacohneg_complete$mod[i] %in% x))]
  if(length(temp) == 1){
    copacohneg_complete$copa_unique[i] <- temp
  } else if(length(temp) == 2){
    copacohneg_complete$copa_unique[i] <- paste(temp, collapse = " ")
  } else {
    copacohneg_complete$copa_unique[i] <- "multiple"
  }
}
copacohneg_complete <- copacohneg_complete |>
  mutate(max_AD = lapply(mod, \(modind){
    mean1 <- sea_means[[1]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
      return(names(mean1)[which.max(mean1)])
    }) |> unlist()) |>
  mutate(max_con = lapply(mod, \(modind){
    mean1 <- sea_means[[2]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
      return(names(mean1)[which.max(mean1)])
    }) |> unlist()) |>
  mutate(dist_all = lapply(mod, \(modind){
    seadist$all[which(these_mods %in% modind)]
    }) |> unlist()) |>
  mutate(dist_per = lapply(mod, \(modind){ # Add info about subclass with highest percentage diff
    tempdf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[modind]], ] |> colMeans() |> unlist()
    return(names(tempdf)[which.min(tempdf)])
    }) |> unlist()) |>
  dplyr::filter((copa_unique == max_AD | copa_unique == max_con) & subclass == copa_unique) |>
  arrange(desc(dist_all))


# Plot percentage differences
sub_diff_mod <- lapply(thesead2, \(mod){
  moddf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[mod]], ] |>
    pivot_longer(everything(), names_to = "subclass", values_to = "pcnt_diff") |>
    mutate("modno" = mod) 
  return(moddf)
}) |> do.call(what = "rbind") |>
  group_by(subclass, modno) |>
  summarise("meanpd" = mean(pcnt_diff)) |>
  mutate("highlight" = F)
for(i in 1:nrow(sub_diff_mod)){
  if(sub_diff_mod$subclass[i] %in% copacohpos_complete$subclass[copacohpos_complete$mod == sub_diff_mod$modno[i]]){
    sub_diff_mod$highlight[i] <- T
  }
}  

sub_diff_mod1 <- lapply(thesead2, \(mod){
  moddf <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], ] |>
    pivot_longer(everything(), names_to = "subclass", values_to = "diff") |>
    mutate("modno" = mod) 
  return(moddf)
}) |> do.call(what = "rbind") |>
  group_by(subclass, modno) |>
  summarise("meanpd" = mean(diff)) |>
  mutate("highlight" = F)
for(i in 1:nrow(sub_diff_mod1)){
  if(sub_diff_mod1$subclass[i] %in% copacohpos_complete$subclass[copacohpos_complete$mod == sub_diff_mod1$modno[i]]){
    sub_diff_mod1$highlight[i] <- T
  }
}  

# Plot neg modules
theseadneg <- unique(copacohneg_complete[,1]) |> as.numeric()
# Which of the coherent mods align with copa mods unique for 1 subclass
eucnegcommon <- unlist(euc_dist$neg[[3]])
eucnegcommon <- unique(eucnegcommon[duplicated(eucnegcommon)])
eucnegu <- lapply(euc_dist$neg[[3]], \(x){
  return(x[!x %in% eucnegcommon])
})
#sum(theseadneg %in% unlist(eucnegu))
# 11 of the 20 coherent modules are significantly unique for a single celltype
theseadneg2 <- theseadneg[theseadneg %in% unlist(eucnegu)]

# Plot percentage differences (neg)
sub_diff_modneg <- lapply(theseadneg2, \(mod){
  moddf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[mod]], ] |>
    pivot_longer(everything(), names_to = "subclass", values_to = "pcnt_diff") |>
    mutate("modno" = mod) 
  return(moddf)
}) |> do.call(what = "rbind") |>
  group_by(subclass, modno) |>
  summarise("meanpd" = mean(pcnt_diff)) |>
  mutate("highlight" = F)
for(i in 1:nrow(sub_diff_modneg)){
  if(sub_diff_modneg$subclass[i] %in% copacohneg_complete$subclass[copacohneg_complete$mod == sub_diff_modneg$modno[i]]){
    sub_diff_modneg$highlight[i] <- T
  }
}  

sub_diff_modneg1 <- lapply(theseadneg2, \(mod){
  moddf <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], ] |>
    pivot_longer(everything(), names_to = "subclass", values_to = "diff") |>
    mutate("modno" = mod) 
  return(moddf)
}) |> do.call(what = "rbind") |>
  group_by(subclass, modno) |>
  summarise("meanpd" = mean(diff)) |>
  mutate("highlight" = F)
for(i in 1:nrow(sub_diff_modneg1)){
  if(sub_diff_modneg1$subclass[i] %in% copacohneg_complete$subclass[copacohneg_complete$mod == sub_diff_modneg1$modno[i]]){
    sub_diff_modneg1$highlight[i] <- T
  }
}  


### Plot modules
# source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_plotting_fxn.R"))
# plot_mods(ind = subdiffmodall$modno,
#           save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels"),
#           file_name = "mod_list.png")
plotmods <- function(subdiffallabs, 
                     subdiffallper, 
                     adconmeans,
                     inds = unique(subdiffallper$modno),
                     save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels"),
                     file_name = "mod_list_all.png"
                     ){
  save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
  module_output_dir <- save_dir1
  all_plots = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_bc.qs")), 
                  qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs")))
  # for indexing and file names
  datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
  if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
  mods <- tapply(datkme[,2], datkme[,3], list)
  modulelengths <- unlist(lapply(mods,length))
  filter_under <- 3
  these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
  cols <- RColorBrewer::brewer.pal(6, "Paired")
  cols2 <- RColorBrewer::brewer.pal(10, "Spectral")

  out_plot <- lapply(seq_along(inds), \(i){
    j <- which(these_mods == inds[i])
    newtheme <- theme_light() + theme(axis.text.x = element_text(size = 6), 
                                    axis.text.y = element_text(size = 6), 
                                    axis.title.x = element_text(size = 6), 
                                    axis.title.y = element_text(size = 6), 
                                    legend.text = element_text(size = 6),
                                    plot.subtitle = element_blank(),
                                    legend.key.size = unit(0.2, "cm"), 
                                    legend.title = element_blank())
    plots <- lapply(all_plots, \(x) x[[j]] + newtheme)

    plots[[1]]$layers[[1]] <- NULL
    plots[[1]] <- plots[[1]] + 
        geom_line(linewidth = 0.2) +
        labs(x = "Sample") + 
        theme(axis.title.x = element_text(size = 6), 
              axis.title.y = element_text(size = 4),
              legend.text = element_text(size = 7), 
              legend.box.margin = margin(0, 0, 0, -10),
              plot.title = element_text(size = 6, hjust = 0.5, margin = margin(0, 0, -0.1, 0))) +
        scale_color_manual(values = cols2) +
        ggtitle(paste0("Module ", inds[i]))
        

    plots[[2]] <- all_plots[[2]][[j]]$data |> 
        dplyr::filter(pval>0) |>
        dplyr::mutate(SetName=factor(SetName, levels=rev(unique(SetName)))) |>
        ggplot(aes(x = SetName, y = pval)) +
            theme_light() +
            geom_bar(stat="identity") +
            theme(axis.text.x = element_text(size=6,hjust=1,vjust = 0.5),
                #axis.title.x = element_text(size=6),
                plot.title = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size=3),
                axis.text.y = element_text(size=4),
                plot.margin = margin(0,0,0,-20)) +
            labs(x="", y=bquote(-log[10](p-val))) #+
            #scale_x_discrete(limits = rev(levels(SetName))) +
            #geom_hline(yintercept = gsea_cutFDR, color = "red")
    plots[[2]]$layers[[2]] <- all_plots[[2]][[j]]$layers[[2]]
    plots[[2]] <- plots[[2]] + coord_flip()
    #p <- plot_grid(plotlist = plots, ncol = 4, nrow = 1, align = "h", axis = "bt", rel_widths=c(0.5, 0.85, 0.75, 1, 0.8))
    #return(p)

    plots[[3]] <- adconmeans |> 
      dplyr::filter(mod == inds[i]) |>
      ggplot(aes(x = subclass, y = mean, fill = type)) + 
        theme_classic() + 
        geom_boxplot(outlier.size = 0.1) +
        theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=6),
              axis.text.y = element_text(size=6),
              legend.position = "none",
              plot.title = element_blank()) +
        labs(y = "mean expression")


    plots[[4]] <- subdiffallabs |> 
      dplyr::filter(modno == inds[i]) |>
      ggplot(aes(x = subclass, y = meanpd, fill = highlight)) +
        theme_classic() + 
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=6),
              axis.text.y = element_text(size=6),
              #legend.text = element_text(size = 6),
              #legend.key.size = unit(0.2, "cm"),
              legend.position = "none",
              plot.title = element_blank()) + 
        labs(y = "mean change in expr")
    
    plots[[5]] <- subdiffallper |> 
      dplyr::filter(modno == inds[i]) |>
      ggplot(aes(x = subclass, y = meanpd, fill = highlight)) +
        theme_classic() + 
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=6),
              axis.text.y = element_text(size=6),
              #legend.text = element_text(size = 6),
              #legend.key.size = unit(0.2, "cm"),
              legend.position = "none",
              plot.title = element_blank()) + 
        labs(y = "mean % change in expr")

    return(plots)
  })

  out_plot2 <- lapply(c(1,2,3,4,5), \(x){
      plist <- lapply(out_plot, \(y) y[[x]])
      if(x != 2){
        return(plot_grid(plotlist = plist, nrow = length(plist), align = "v", axis = "rl"))
      } else {
        return(plot_grid(plotlist = plist, nrow = length(plist)))
      }
  })

  outall2 <- plot_grid(plotlist = out_plot2, ncol = 5, align = "h", axis = "bt", rel_widths = c(0.6, 0.4, 0.8, 0.8, 0.8))
  ggsave(outall2, file = file.path(save_path, file_name),  height = length(inds) * 1, width = 12, bg = "white", limitsize = F)
}

subdiffallper <- rbind(sub_diff_mod, sub_diff_modneg) |>
  mutate(meanpd = meanpd * 100)
subdiffallabs <- rbind(sub_diff_mod1, sub_diff_modneg1)
adconmeans <- lapply(unique(subdiffallper[[2]]), \(mod){
  mapply(\(sea, name){
    sea[rownames(sea) %in% mod_bc[[mod]], ] |>
      pivot_longer(everything(), names_to = "subclass", values_to = "mean") |>
      mutate("type" = name)
  }, sea_means[1:2], c("AD", "Con"), SIMPLIFY = F) |>
    do.call(what = "rbind") |>
    mutate(mod = mod)
}) |> do.call(what = "rbind")

plotmods(subdiffallabs = subdiffallabs,
         subdiffallper = subdiffallper, 
         adconmeans = adconmeans,
         inds = unique(subdiffallper$modno),
         save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels"),
         file_name = "mod_list_all.png")

inds_filt <- c(unique(copacohpos_complete$mod, copacohneg_complete$mod))

plotmods(subdiffallabs = subdiffallabs,
         subdiffallper = subdiffallper, 
         adconmeans = adconmeans,
         inds = inds_filt, # all copacohneg mods filtered out
         save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels"),
         file_name = "mod_list_all_filtered_sorted.png")

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_expression_consistency_fxns.R"))
find_cons_mods(save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/"))


##################################################
# Let's add early/late ad modules to this calculus
##################################################

euc_dist_el <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/euclidean_distances/euclidean_sigmods_all.qs")),
                    "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/euclidean_distances/euclidean_sigmods_positive.qs")),
                    "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/euclidean_distances/euclidean_sigmods_negative.qs")))
sea_means_el <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/sn_summary_tables/allmlist_log.qs"))
sub_diff_el <- sea_means_el[[1]] - sea_means_el[[2]]

ccpos <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff_el[rownames(sub_diff_el) %in% mod_bc[[mod]], colnames(sub_diff_el) == subclass_name]
    return(sum(temp > 0) / length(temp))
  }) |> unlist() 
  if(length(outvec) > 0){
    return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
  } else {
    return(data.frame("subclass" = subclass_name, "pcnt" = NA))
  }
}, euc_dist_el$pos[[3]], names(euc_dist_el$pos[[3]]), SIMPLIFY = F)
# How many of these overlap with the conVsAD modules?
t1 <- lapply(ccpos, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # earlyvslate
t2 <- lapply(copacohpos, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # conVsAD
unique(t1)[which(unique(t1) %in% unique(t2))]
# 149 716
unique(t2)[which(unique(t2) %in% unique(t1))]
# 149 716

# Which are specific to a subclass?
ccpos_com <- lapply(ccpos, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
#tapply(ccpos_com[ ,2], ccpos_com[ ,1], list)
ccpos_u <- tapply(ccpos_com[ ,2], ccpos_com[ ,1], list)
ccpos_u <- ccpos_u[(lapply(ccpos_u, length) |> unlist()) == 1] |> unlist()
    #      25         295         301         441         716         744 
    # "L6 CT"     "L6 CT"     "L6 CT"      "Sncg"     "L6 CT"       "Vip" 
    #      87         909 
    #  "VLMC" "Astrocyte" 

# Negative (higher in con)
ccneg <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff_el[rownames(sub_diff_el) %in% mod_bc[[mod]], colnames(sub_diff_el) == subclass_name]
    return(sum(temp < 0) / length(temp))
  }) |> unlist() 
  return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
}, euc_dist_el$neg[[3]], names(euc_dist_el$neg[[3]]), SIMPLIFY = F) 
lapply(ccneg, \(x) sum(x$pcnt == 1)) |> unlist()
#       Astrocyte      Chandelier     Endothelial         L2/3 IT           L4 IT 
#               2               7               1               1              37 
#           L5 ET           L5 IT         L5/6 NP           L6 CT           L6 IT 
#               0               5               1               2               3 
#      L6 IT Car3             L6b           Lamp5      Lamp5 Lhx6   Microglia-PVM 
#               7               2              18               2               0 
# Oligodendrocyte             OPC            Pax6           Pvalb            Sncg 
#              15               3               1               1               1 
#             Sst       Sst Chodl             Vip            VLMC 
#               0               0               1               0 
# How many of these overlap with the conVsAD modules?
t1 <- lapply(ccneg, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # earlyvslate
t2 <- lapply(copacohneg, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # conVsAD
unique(t1)[which(unique(t1) %in% unique(t2))]
# quite a few
unique(t2)[which(unique(t2) %in% unique(t1))]
# quite a few
# How many of these modules are specific to a single subclass?
ccneg_com <- lapply(ccneg, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
ccneg_u <- tapply(ccneg_com[ ,2], ccneg_com[ ,1], list)
ccneg_u <- ccneg_u[(lapply(ccneg_u, length) |> unlist()) == 1] |> unlist()
#                10               101              1035              1049                                                                                         
#           "L4 IT"           "L4 IT"           "L6 CT"       "Astrocyte" 
#              1086                11                13               131 
#           "L4 IT"           "L4 IT"           "L4 IT"           "L4 IT" 
#                14               140               141               147 
#           "L4 IT"           "Lamp5"           "L4 IT"       "Astrocyte" 
#               157                17               173               184 
#           "L4 IT"           "Lamp5" "Oligodendrocyte"           "Lamp5" 
#               192                21               227                25 
# "Oligodendrocyte"           "L4 IT"           "L4 IT" "Oligodendrocyte" 
#               260               272                28               316 
# "Oligodendrocyte"           "L4 IT"           "L4 IT"           "L4 IT" 
#               329               334                35               442 
#      "L6 IT Car3"           "L4 IT"           "L4 IT"         "L2/3 IT" 
#               452               453                48               489 
#           "Lamp5"            "Sncg"           "L4 IT"           "Lamp5" 
#               493               517               518               526 
#      "L6 IT Car3"           "Lamp5"      "L6 IT Car3" "Oligodendrocyte" 
#                58               614               623               639 
#           "Lamp5"           "L6 IT" "Oligodendrocyte"      "Chandelier" 
#                66               678               716               721 
#           "L4 IT"           "L4 IT" "Oligodendrocyte"      "L6 IT Car3" 
#               754                82                88               893 
# "Oligodendrocyte"           "L4 IT"      "Chandelier"     "Endothelial" 
#                90               910               968               984 
#           "L4 IT"           "L6 CT"           "L4 IT"         "L5/6 NP" 


# Let's add con/early
euc_dist_ce <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/euclidean_distances/euclidean_sigmods_all.qs")),
                    "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/euclidean_distances/euclidean_sigmods_positive.qs")),
                    "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/euclidean_distances/euclidean_sigmods_negative.qs")))
sea_means_ce <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
sub_diff_ce <- sea_means_ce[[1]] - sea_means_ce[[2]]

ccposce <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff_ce[rownames(sub_diff_ce) %in% mod_bc[[mod]], colnames(sub_diff_ce) == subclass_name]
    return(sum(temp > 0) / length(temp))
  }) |> unlist() 
  if(length(outvec) > 0){
    return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
  } else {
    return(data.frame("subclass" = subclass_name, "pcnt" = NA))
  }
}, euc_dist_ce$pos[[3]], names(euc_dist_ce$pos[[3]]), SIMPLIFY = F)
# # How many of these overlap with the conVsAD modules?
# t1 <- lapply(ccposce, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # convsearly
# t2 <- lapply(copacohpos, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # conVsAD
# unique(t1)[which(unique(t1) %in% unique(t2))]
# # 948 149
# unique(t2)[which(unique(t2) %in% unique(t1))]
# # 149 948

# Which are specific to a subclass?
ccposce_com <- lapply(ccposce, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
#tapply(ccpos_com[ ,2], ccpos_com[ ,1], list)
ccposce_u <- tapply(ccposce_com[ ,2], ccposce_com[ ,1], list)
ccposce_u <- ccposce_u[(lapply(ccposce_u, length) |> unlist()) == 1] |> unlist()
    #      25         295         301         441         716         744 
    # "L6 CT"     "L6 CT"     "L6 CT"      "Sncg"     "L6 CT"       "Vip" 
    #      87         909 
    #  "VLMC" "Astrocyte" 

# Negative (higher in con)
ccnegce <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff_ce[rownames(sub_diff_ce) %in% mod_bc[[mod]], colnames(sub_diff_ce) == subclass_name]
    return(sum(temp < 0) / length(temp))
  }) |> unlist() 
  return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
}, euc_dist_ce$neg[[3]], names(euc_dist_ce$neg[[3]]), SIMPLIFY = F) 
lapply(ccnegce, \(x) sum(x$pcnt == 1)) |> unlist()
#       Astrocyte      Chandelier     Endothelial         L2/3 IT           L4 IT 
#              19               4              12               5               5 
#           L5 ET           L5 IT         L5/6 NP           L6 CT           L6 IT 
#               3              17               1               0              11 
#      L6 IT Car3             L6b           Lamp5      Lamp5 Lhx6   Microglia-PVM 
#              18              36               9               6               7 
# Oligodendrocyte             OPC            Pax6           Pvalb            Sncg 
#              23               6               0               7               3 
#             Sst       Sst Chodl             Vip            VLMC 
#              10               0               4               0
# How many of these overlap with the conVsAD modules?
t1 <- lapply(ccnegce, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # convsearly
t2 <- lapply(copacohneg, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # conVsAD
unique(t1)[which(unique(t1) %in% unique(t2))]
# many
unique(t2)[which(unique(t2) %in% unique(t1))]
# many
# How many of these modules are specific to a single subclass?
ccnegce_com <- lapply(ccnegce, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
ccnegce_u <- tapply(ccnegce_com[ ,2], ccnegce_com[ ,1], list)
ccnegce_u <- ccnegce_u[(lapply(ccnegce_u, length) |> unlist()) == 1] |> unlist()

## Flowchart
# Find modules =>
# ...that are significant (COPA) =>
# ...and where all mod genes are moving in same direction for any given celltype =>
# ...and where the significant subclass (COPA) is the same as the subclass where all genes are moving in same direction =>
# ...and where that subclass is also the subclass with highest average expression for that module
 
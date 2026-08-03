# Table structure:
#              Type:    |      Ctrl vs all AD  |  Control vs Early                |  Early vs Late
#                           Up in AD | down in AD:   Up in Early | down in Early:        Up in Late | down in Late:
# Dataset:
# SEAAD2024                     10          2            36               2                0              14
# Morabito [other AD dataset]   3          5
# Rosmap AD [bulk AD mods]      14         19   

library(data.table)
library(qs)
library(tidyverse)
library(gt)
library(gtExtras)
library(webshot2) # saving gt tables as image
library(showtext)
showtext_auto()

#############################
# Load and count objects (AD)
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
      outlist <- lapply(df, \(a) a$mod)
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

############
# Load modules
############
# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Load expr (ROSMAP bulk AD samples only)
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_AD_only_SampleNetworks/1_10-52-54/ROSMAP_samp_filt_AD_only_1_248_ComBat.csv"), data.table=F)
megaset_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
expr <- expr[match(megaset_expr[,1], expr[,1]),]
expr <- data.frame("ensembl_id"=expr[,1],"Gene"= megaset_expr[,2], expr[,3:ncol(expr)])
expr <- expr[!is.na(expr[,1]),]
bulk_expr_ros <- expr
bulkt_ros <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
# Load mods (ROSMAP bulk AD samples only)
mod_seed_ros <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024/modules/unmerged_modules.qs"))
mod_bc_ros <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_ros_lengths <- lapply(mod_bc_ros, length) |> unlist()

###############
# Fisher's exact for significant bulk vs ROSMAP modules
#############
fisherTest=function(set,mod,all){
  total.shared=length(intersect(all,set))
  shared.in.mod=length(intersect(mod,set))
  shared.out.mod=total.shared-shared.in.mod
  in.mod.not.shared=length(mod)-shared.in.mod
  out.mod.not.shared=length(all)-length(mod)-shared.out.mod
  fisher.test(matrix(c(shared.in.mod,in.mod.not.shared,shared.out.mod,out.mod.not.shared),ncol=2),alternative="greater")$p.val
}

all_genes <- intersect(bulk_expr$Gene, bulk_expr_ros$Gene)

# bulk mods vs rosmap mods
testlist <- list()
for(i in seq_along(mod_bc)){
  modvec <- c()
  for(j in seq_along(mod_bc_ros)){
    modvec[j] <- fisherTest(mod_bc_ros[[j]], mod_bc[[i]], all_genes)
  }
  testlist[[i]] <- modvec
  cat(i, " ")
}
testdf <- do.call(cbind, testlist)
fwrite(testdf, file = "~/test/test.csv")
testdf <- fread(data.table = F, file = "~/test/test.csv")
# rows: rosmapmods
# cols: bulk mods
#testdfcolranks <- apply(testdf, 2, rank)
 
# bulk mods vs bulk mods
testlist2 <- list()
for(i in seq_along(mod_bc)){
  modvec <- c()
  for(j in seq_along(mod_bc)){
    modvec[j] <- fisherTest(mod_bc[[j]], mod_bc[[i]], all_genes)
  }
  testlist2[[i]] <- modvec
  cat(i, " ")
}
testdf2 <- do.call(cbind, testlist2)
fwrite(testdf2, file = "~/test/test2.csv")
testdf2 <- fread(data.table = F, file = "~/test/test2.csv")

# Which rosmap mods correspond to bulk megaset mods
#lapply(pathList$bulk_megaset$conVAll$pos, \(x) which.min(testdf[, as.numeric(x)])) |> unlist() # which rosmap mod has lowest pval
#lapply(pathList$bulk_megaset$conVAll$pos, \(x) min(testdf[, as.numeric(x)])) |> unlist() # what is lowest pval
#lapply(pathList$bulk_megaset$conVAll$pos, \(x) testdfcolranks[as.numeric(pathList$rosmapAD$conVsAllAD$pos), as.numeric(x)]) # what are ranks of significant rosmap mods for bulk mods
library(ComplexHeatmap)
pdf("~/test/test.pdf")
for(i in 1:4){
  bulk1 <- pathList$bulk_megaset[[i]] |> unlist()
  ros1 <- pathList$rosmapAD[[i]] |> unlist()
  plotdf <- testdf[as.numeric(ros1), as.numeric(bulk1)]
  draw(Heatmap(plotdf, row_title= "Significant mods (ROSMAP SEAAD2024)", column_title = "Significant mods (Bulk SEAAD2024)"))
}
dev.off()

# Which bulk mods (seaad2024) correspond to bulk mods (mit)
pdf("~/test/test2.pdf")
for(i in 1:4){
  bulk1 <- pathList$bulk_megaset[[i]] |> unlist()
  ros1 <- pathList$bulk_MIT[[i]] |> unlist()
  plotdf <- testdf2[as.numeric(ros1), as.numeric(bulk1)]
  draw(Heatmap(plotdf, row_title= "Significant mods (MIT)", column_title = "Significant mods (SEAAD2024)"))
}
dev.off()

# Which modules are conserved?
# - In bulk mods vs rosmap on SEAAD2024?
inros <- list()
inmit <- list()
inames <- c("Con vs all AD", "Con vs early AD", "Con vs late AD", "Early vs late AD")
for(i in 1:4){
  bulk1 <- pathList$bulk_megaset[[i]] |> unlist()
  mit1 <- pathList$bulk_MIT[[i]] |> unlist()
  ros1 <- pathList$rosmapAD[[i]] |> unlist()
  plotdf <- testdf[as.numeric(ros1), as.numeric(bulk1)]
  plotdf2 <- testdf2[as.numeric(mit1), as.numeric(bulk1)]

  conserved_mods <- apply(plotdf, 2, \(x) sum(x < 0.05))
  conserved_mods2 <- apply(plotdf2, 2, \(x) sum(x < 0.05))
  conserved_mods3 <- apply(plotdf2, 1, \(x) sum(x < 0.05))

  inros[[i]] <- data.frame("comp" = inames[i],
                           "which" = c("significant using bulk megaset mods", "significant using ROSMAP mods\nor bulk megaset mods"),
                           "val" = c(sum(conserved_mods == 0), sum(conserved_mods > 0)))
  inmit[[i]] <- data.frame("comp" = inames[i],
                           "which" = c("significant when projected on SEAAD2024",
                                       "significant when projected on MIT_Multiome", 
                                       "significant for both"),
                           "val" = c(sum(conserved_mods2 == 0), sum(conserved_mods3 == 0), sum(conserved_mods2 > 0)))
}

inrosdf <- do.call(rbind, inros)
p <- ggplot(inrosdf, aes(x = comp, y = val, fill = which)) +
  theme_minimal() + 
  geom_bar(stat = "identity", position = "stack", alpha = 0.5) +
  labs(x = "", y = "# of significant modules") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        legend.title = element_blank())
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/fisher_overlap_bulkMegasetSEA_vs_rosmapSEA.svg"), bg = "white", width = 8)

inmitdf <- do.call(rbind, inmit)
p <- inmitdf |>
  filter(which == "significant for both") |>
  ggplot(aes(x = comp, y = val, fill = which)) +
  theme_minimal() + 
  geom_bar(stat = "identity", position = "stack", alpha = 0.5) +
  labs(x = "", y = "# of shared significant modules") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        legend.title = element_blank(),
        legend.position = "none") +
        paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/fisher_overlap_bulkMegasetSEA_vs_bulkMegasetMIT.svg"), bg = "white", width = 5)

## Plot celltypes
pathList_ct <- lapply(paths[1], \(d){
  compname <- list.files(d)
  modlist <- lapply(compname, \(x){
    path <- file.path(d, x, "modlist_filtered.qs")
    if(file.exists(path)){
      df <- qread(file.path(d, x, "modlist_filtered.qs"))
      outlist <- mapply(\(a,b){
        out <- a$subclass |> table() |> data.frame()
        out$name <- b
        return(out)
      }, df, c("pos", "neg"), SIMPLIFY = F) |> do.call(what = "rbind")
    } else {
      outlist <- list()
    }
    outlist$comp <- x
    return(outlist)
  })
  return(modlist |> do.call(what = "rbind"))
})[[1]]


cols <- RColorBrewer::brewer.pal(length(unique(pathList_ct$Var1)), "Set3")
names(cols) = unique(pathList_ct$Var1)
colsdf <- data.frame("Var1" = unique(pathList_ct$Var1), "cols" = cols)

comps <- unique(pathList_ct$comp)
nameList <- list(c("higher in AD", "higher in Con"),
                 c("higher in Con", "higher in early AD"),
                 c("higher in Con", "higher in late AD"),
                 c("higher in late AD", "higher in early AD"))
for(i in 1:4){
  plotdf <- pathList_ct |> filter(comp == comps[i])
  plotdf$name2 <- ifelse(plotdf$name == "pos", nameList[[i]][1], nameList[[i]][2])
  plotdf$Var1 <- factor(plotdf$Var1, levels = unique(plotdf$Var1))
  if(i == 4){
    plotdf$name2 <- factor(plotdf$name2, levels = rev(unique(plotdf$name2)))
  } else {
    plotdf$name2 <- factor(plotdf$name2, levels = unique(plotdf$name2))
  }
  plotdf <- left_join(plotdf, colsdf)

  p <- ggplot(plotdf, aes(x = name2, y = Freq, fill = Var1)) +
    theme_minimal() + 
    geom_bar(stat = "identity", position = "stack", alpha = 0.9) +
    labs(x = "", y = "# of DECoPA modules") +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          plot.margin = margin(1, 1, 1, 3, "cm"), 
          legend.title = element_blank()) +
    scale_fill_manual(values = cols)
  ggsave(p, file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/stackbarplot_celltypes_bulkmegasetOnSea_"), comps[i], ".svg"), width = 6, bg = "white")
}

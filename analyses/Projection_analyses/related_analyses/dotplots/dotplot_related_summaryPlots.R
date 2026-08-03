library(tidyverse)
library(data.table)
library(qs)
library(ggpubr)
library(UpSetR)
library(ComplexHeatmap)
library(showtext)
showtext_auto()

# First, load all data:
dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/")

title_vec <- c("Con vs all AD (DFC)",
             "Early AD vs Con (DFC)",
             "Late vs Early AD (DFC)",
             "Con vs all AD (MTG)",
             "Early AD vs Con (MTG)",
             "Late vs Early AD (MTG)",
             "Con vs all AD (ROSMAP, DFC)",
             "Early AD vs Con (ROSMAP, DFC)",
             "Late vs Early AD (ROSMAP, DFC)",
             "Con vs all AD (ROSMAP, MTG)",
             "Early AD vs Con (ROSMAP, MTG)",
             "Late vs Early AD (ROSMAP, MTG)",
             "APOE 4/4 vs 3/3 (DFC)",
             "SCZ vs con (DFC)")

dat_vec <- c("gabitto_vs_liu_DFC_AllADVsCon",    
             "gabitto_vs_liu_DFC_earlyVsCon",        
             "gabitto_vs_liu_DFC_lateVsEarly",        
             "gabitto_vs_liu_MTG_AllADVsCon",         
             "gabitto_vs_liu_MTG_earlyVsCon",   
             "gabitto_vs_liu_MTG_lateVsEarly",        
             "gabitto_vs_liu_DFC_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_DFC_earlyVsCon_ROSMAP", 
             "gabitto_vs_liu_DFC_lateVsEarly_ROSMAP",
             "gabitto_vs_liu_MTG_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_MTG_earlyVsCon_ROSMAP", 
             "gabitto_vs_liu_MTG_lateVsEarly_ROSMAP",
             "gabitto_vs_liu_APOE_DFC",
             "brainSCOPE_CMC_vs_SZBD") 

label_vec <- c("Shared in Gabitto and Liu DFC (con vs all AD)",    
               "Shared in Gabitto and Liu DFC (early AD vs con)",        
               "Shared in Gabitto and Liu DFC (late vs early AD)",        
               "Shared in Gabitto and Liu MTG (con vs all AD)",         
               "Shared in Gabitto and Liu MTG (early AD vs con)",   
               "Shared in Gabitto and Liu MTG (late vs early AD)",        
               "Shared in Gabitto and Liu DFC (con vs all AD, ROSMAP mods)",    
               "Shared in Gabitto and Liu DFC (early AD vs con, ROSMAP mods)", 
               "Shared in Gabitto and Liu DFC (late vs early AD, ROSMAP mods)",
               "Shared in Gabitto and Liu MTG (con vs all AD, ROSMAP mods)",    
               "Shared in Gabitto and Liu MTG (early AD vs con, ROSMAP mods)", 
               "Shared in Gabitto and Liu MTG (late vs early AD, ROSMAP mods)",
               "Shared in Gabitto and Liu DFC (APOE 4/4 vs 3/3)",
               "Shared in brainSCOPE CMC and SZBD DFC (SCZ vs con)") 

label_vec_abrv <- c("DFC (con vs all AD)",    
                    "DFC (early AD vs con)",        
                    "DFC (late vs early AD)",        
                    "MTG (con vs all AD)",         
                    "MTG (early AD vs con)",   
                    "MTG (late vs early AD)",        
                    "DFC (con vs all AD, ROSMAP mods)",    
                    "DFC (early AD vs con, ROSMAP mods)", 
                    "DFC (late vs early AD, ROSMAP mods)",
                    "MTG (con vs all AD, ROSMAP mods)",    
                    "MTG (early AD vs con, ROSMAP mods)", 
                    "MTG (late vs early AD, ROSMAP mods)",
                    "DFC (APOE 4/4 vs 3/3)",
                    "DFC (SCZ vs con)") 

dat_list <- mapply(\(x, y){
  out <- fread(file.path(dir1, x, "shared_output_table.csv"), data.table = F)
  out <- out[!duplicated(paste0(out$mod, out$Celltype)), ] |>
    dplyr::mutate(Direction = factor(Direction, levels = c(-1, 1)))
  l_out <- tapply(out$mod, out$Direction, list)
  return(l_out)
}, dat_vec, title_vec, SIMPLIFY = F)
names(dat_list) <- label_vec

dat_list_neg <- lapply(dat_list, \(x) x[[1]])
dat_list_pos <- lapply(dat_list, \(x) x[[2]])

# 1. Find the modules that overlap between comparisons and create...
#  - ...Upset plot
#  - ...Heatmap
#  - ...Node-link network
# - Add hypergeom p-values
p <- upset(fromList(dat_list_neg), 
      sets = names(dat_list_neg), 
      order.by = "freq",
      mainbar.y.label = "# of overlapping modules", 
      sets.x.label = "# of modules (Lower in severe cases)")

pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/upset_lowerInMoreSevere.pdf"), width = 16)
p
dev.off()

p <- upset(fromList(dat_list_pos), 
      order.by = "freq",
      sets = names(dat_list_pos), 
      mainbar.y.label = "# of overlapping modules", 
      sets.x.label = "# of modules (Higher in severe cases)")
pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/upset_higherInMoreSevere.pdf"), width = 16)
p
dev.off()

# Heatmap (of Jaccard index)

make_jac_mat <- function(set_list){
  num_sets <- length(set_list)
  jaccard_matrix <- matrix(NA, nrow = num_sets, ncol = num_sets)

  for (i in 1:num_sets) {
    for (j in 1:num_sets) {
      intersection_size <- length(intersect(set_list[[i]], set_list[[j]]))
      union_size <- length(union(set_list[[i]], set_list[[j]]))
      jaccard_matrix[i, j] <- intersection_size / union_size
    }
  }
  rownames(jaccard_matrix) <- label_vec_abrv
  colnames(jaccard_matrix) <- label_vec_abrv
  return(jaccard_matrix)
}

jac_neg <- make_jac_mat(dat_list_neg)
jac_pos <- make_jac_mat(dat_list_pos)

dim1 <- 8
m1 <- 3

pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/heatmap_lowerInMoreSevere.pdf"), width = dim1, height = dim1/1.3)
p <- Heatmap(jac_neg,
             show_column_names = F,
             heatmap_legend_param = list(title = "Jaccard index"),
             column_title = "Mods shared between Gabitto 2024 and Liu 2025\n(lower in more severe samples)")
draw(p, 
     padding = unit(c(m1, m1, m1, 2 * m1), "cm"),
     heatmap_legend_side = "bottom")
dev.off()

pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/heatmap_higherInMoreSevere.pdf"), width = dim1, height = dim1/1.3)
p <- Heatmap(jac_pos,
             show_column_names = F,
             heatmap_legend_param = list(title = "Jaccard index"),
             column_title = "Mods shared between Gabitto 2024 and Liu 2025\n(higher in more severe samples)")
draw(p, 
     padding = unit(c(m1, m1, m1, 2 * m1), "cm"),
     heatmap_legend_side = "bottom")
dev.off()

# Network graph

library(igraph)
library(ggraph)

# From distance matrix (convert to similarity: 1-dist)
g <- graph_from_adjacency_matrix(jac_neg, mode="undirected", weighted=TRUE)

# Filter edges for better visualization (e.g., keep only strong similarities)
# g_filtered <- delete_edges(g, E(g)[weight < threshold])

p <- ggraph(g, layout = 'kk') + # Or other layouts like 'fr', 'circle'
  geom_edge_fan(aes(width = weight, alpha = weight), color = "gray50") + # Weight/alpha by Jaccard
  geom_node_point(size = 5, color = "steelblue") +
  geom_node_label(aes(label = name), repel = TRUE, nudge_x = 0, min.segment.length = 0.01, force_pull = 2, alpha = 0.8) +
  theme_graph() +
  labs(edge_width = "Jaccard Similarity", edge_alpha = "Jaccard Similarity")

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/network_lowerInMoreSevere.pdf"))

g <- graph_from_adjacency_matrix(jac_pos, mode="undirected", weighted=TRUE)

# Filter edges for better visualization (e.g., keep only strong similarities)
# g_filtered <- delete_edges(g, E(g)[weight < threshold])

p <- ggraph(g, layout = 'kk') + # Or other layouts like 'fr', 'circle'
  geom_edge_fan(aes(width = weight, alpha = weight), color = "gray50") + # Weight/alpha by Jaccard
  geom_node_point(size = 5, color = "steelblue") +
  geom_node_label(aes(label = name), repel = TRUE, nudge_x = 0, min.segment.length = 0.01, force_pull = 2, alpha = 0.8) +
  theme_graph() +
  labs(edge_width = "Jaccard Similarity", edge_alpha = "Jaccard Similarity")

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/network_higherInMoreSevere.pdf"))

#############
# Gene overlap analysis
############

# What is hypergeometric p-value of overlap of all significant genes for a given subclass/direction?
# Example: what is p-val for overlap of sig genes for ROSMAP astrocytes (pos) vs bulk megaset astrocytes (pos)?

# hypergeom function:
fisherTest_modoverlap <- function(dat1, # ex. SEAAD2024
                                  dat2, # ex. MIT
                                  shared.in.mod, # intersection between SEAAD2024 and MIT
                                  all){ # Fisher's test function
  #total.shared = length(intersect(all,dat1))
  #shared.in.mod = length(intersect(dat2,dat1))
  shared.out.mod = length(dat1) - shared.in.mod
  in.mod.not.shared = length(dat2) - shared.in.mod
  out.mod.not.shared = length(all) - length(dat2) - shared.out.mod
  fisher.test(matrix(c(shared.in.mod,
                      in.mod.not.shared,
                      shared.out.mod,
                      out.mod.not.shared), ncol=2),alternative="greater")$p.val
}

# Start with bulk megaset vs rosmap
dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/")

class_info <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |>
  arrange(Subclass) |>
  mutate(Subclass_fixed = factor(c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
                                   "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip"),
                                 levels = c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
                                             "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
                                             "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC"))) |>
  arrange(Subclass_fixed) |>
  select(Subclass, Class)

dat_vec <- c("gabitto_vs_liu_DFC_AllADVsCon",    
             "gabitto_vs_liu_DFC_earlyVsCon",        
             "gabitto_vs_liu_DFC_lateVsEarly",        
             "gabitto_vs_liu_MTG_AllADVsCon",         
             "gabitto_vs_liu_MTG_earlyVsCon",   
             "gabitto_vs_liu_MTG_lateVsEarly",        
             "gabitto_vs_liu_DFC_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_DFC_earlyVsCon_ROSMAP", 
             "gabitto_vs_liu_DFC_lateVsEarly_ROSMAP",
             "gabitto_vs_liu_MTG_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_MTG_earlyVsCon_ROSMAP", 
             "gabitto_vs_liu_MTG_lateVsEarly_ROSMAP") 

title_vec <- c("Con vs all AD (DFC)",
             "Early AD vs Con (DFC)",
             "Late vs Early AD (DFC)",
             "Con vs all AD (MTG)",
             "Early AD vs Con (MTG)",
             "Late vs Early AD (MTG)",
             "Con vs all AD (ROSMAP, DFC)",
             "Early AD vs Con (ROSMAP, DFC)",
             "Late vs Early AD (ROSMAP, DFC)",
             "Con vs all AD (ROSMAP, MTG)",
             "Early AD vs Con (ROSMAP, MTG)",
             "Late vs Early AD (ROSMAP, MTG)")

label_vec <- c("Shared in Gabitto and Liu DFC (con vs all AD)",    
               "Shared in Gabitto and Liu DFC (early AD vs con)",        
               "Shared in Gabitto and Liu DFC (late vs early AD)",        
               "Shared in Gabitto and Liu MTG (con vs all AD)",         
               "Shared in Gabitto and Liu MTG (early AD vs con)",   
               "Shared in Gabitto and Liu MTG (late vs early AD)",        
               "Shared in Gabitto and Liu DFC (con vs all AD, ROSMAP mods)",    
               "Shared in Gabitto and Liu DFC (early AD vs con, ROSMAP mods)", 
               "Shared in Gabitto and Liu DFC (late vs early AD, ROSMAP mods)",
               "Shared in Gabitto and Liu MTG (con vs all AD, ROSMAP mods)",    
               "Shared in Gabitto and Liu MTG (early AD vs con, ROSMAP mods)", 
               "Shared in Gabitto and Liu MTG (late vs early AD, ROSMAP mods)") 

dat_list <- mapply(\(x, y){
  out <- fread(file.path(dir1, x, "shared_output_table.csv"), data.table = F)
  return(out)
}, dat_vec, title_vec, SIMPLIFY = F)
names(dat_list) <- label_vec

# Load module data (topmodposbc - all dCoPA analyses done using topmodposbc definitions)
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
filter_under <- 3
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])

comp_names <- title_vec[1:6]

genedf <- lapply(seq_along(comp_names), \(i){
  d1 <- dat_list[[i]]
  d2 <- dat_list[[i + 6]]
  common_cts <- unique(c(d1[,2], d2[,2]))

  outdf <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(common_cts)){
      m1 <- d1 |> filter(Celltype_Gabitto_2024 == common_cts[c],
                        Direction == d)
      g1 <- unique(unlist(mods[m1$mod]))
      m2 <- d2 |> filter(Celltype_Gabitto_2024 == common_cts[c],
                        Direction == d)
      g2 <- unique(unlist(mods[m2$mod]))

      ngene <- 18913 * 24 * 2
      fisher_pval <- fisherTest_modoverlap(g1, 
                                           g2, 
                                           length(unique(intersect(g1, g2))), # intersection 
                                           all = 1:ngene)

      sublist[[c]] <- data.frame("Comparison" = comp_names[i],
                                 "Direction" = d,
                                 "Celltype" = common_cts[c],
                                 "Length1" = length(g1),
                                 "Length2" = length(g2),
                                 "NMods1" = length(unique(m1$mod)),
                                 "NMods2" = length(unique(m2$mod)),
                                 "LengthIntersectGene" = length(unique(intersect(g1, g2))),
                                 "LengthUnionGene" = length(unique(union(g1, g2))), 
                                 "LengthIntersectMod" = length(unique(intersect(m1$mod, m2$mod))),
                                 "Fisher_p_val" = fisher_pval)
    }
    return(sublist |> do.call(what = "rbind"))
  }) |> do.call(what = "rbind") |>
    mutate(log_pval = -log10(Fisher_p_val),
           log_pval = case_match(log_pval, Inf ~ 324, .default = log_pval),
           Direction = factor(Direction, levels = c("-1", "1")),
           Direction = fct_recode(Direction,
                                  "Lower in more severe" = "-1",
                                  "Higher in more severe" = "1"),
           Celltype = factor(Celltype, levels = class_info$Subclass)
           )
  return(outdf)
}) |> do.call(what = "rbind")

p <- ggplot(genedf, aes(x = Celltype, y = log_pval, fill = Direction)) +
  theme_classic() + 
  geom_bar(position = "dodge", stat = "identity") +
  labs(x = "", y = bquote(-log[10]~"p-value")) +
  #geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~Comparison, nrow = 2, ncol = 3, dir = "v", scales = "free_x")

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/gene_overlap_BulkMegasetVsROSMAP_hypergeomPval.pdf"), width = 12, height = 5)

p2 <- ggplot(genedf, aes(x = Celltype, y = LengthUnionGene, fill = Direction)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "", y = "Total # of genes comprising sig mods") +
  #geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~Comparison, nrow = 2, ncol = 3, dir = "v", scales = "free_x")
ggsave(p2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/gene_overlap_BulkMegasetVsROSMAP_totalgenecount.pdf"), width = 12, height = 5)

## Plot overlap of significant genes between subclasses

genedf_overlap <- lapply(seq_along(comp_names), \(i){
  d1 <- dat_list[[i]]
  d2 <- dat_list[[i + 6]]
  common_cts <- unique(c(d1[,2], d2[,2]))

  outdf <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(common_cts)){
      m1 <- d1 |> filter(Celltype_Gabitto_2024 == common_cts[c],
                        Direction == d)
      g1 <- unique(unlist(mods[m1$mod]))
      m2 <- d2 |> filter(Celltype_Gabitto_2024 == common_cts[c],
                        Direction == d)
      g2 <- unique(unlist(mods[m2$mod]))
      
      outvec <- unique(intersect(g1, g2))
      if(length(outvec) > 0){
        sublist[[c]] <- outvec
      } else {
        sublist[[c]] <- "none"
      }
    }
    names(sublist) <- common_cts
    return(sublist)
  }) 
  names(outdf) <- c(-1, 1)
  return(outdf)
})
names(genedf_overlap) <- comp_names

# First create an upset plot for conVAll DFC, lower in severe
p1 <- upset(fromList(genedf_overlap[[1]][[1]]), 
      sets = names(genedf_overlap[[1]][[1]]), 
      order.by = "freq",
      mainbar.y.label = "# of overlapping genes", 
      sets.x.label = "# of genes",
      mb.ratio = c(0.5, 0.5))
p2 <- upset(fromList(genedf_overlap[[1]][[2]]), 
            sets = names(genedf_overlap[[1]][[2]]), 
            order.by = "freq",
            mainbar.y.label = "# of overlapping genes", 
            sets.x.label = "# of genes",
            mb.ratio = c(0.5, 0.5))
p3 <- upset(fromList(genedf_overlap[[2]][[1]]), 
      sets = names(genedf_overlap[[2]][[1]]), 
      order.by = "freq",
      mainbar.y.label = "# of overlapping genes", 
      sets.x.label = "# of genes",
      mb.ratio = c(0.5, 0.5))
p4 <- upset(fromList(genedf_overlap[[2]][[2]]), 
            sets = names(genedf_overlap[[2]][[2]]), 
            order.by = "freq",
            mainbar.y.label = "# of overlapping genes", 
            sets.x.label = "# of genes",
            mb.ratio = c(0.5, 0.5))
p5 <- upset(fromList(genedf_overlap[[3]][[1]]), 
      sets = names(genedf_overlap[[3]][[1]]), 
      order.by = "freq",
      mainbar.y.label = "# of overlapping genes", 
      sets.x.label = "# of genes",
      mb.ratio = c(0.5, 0.5))
p6 <- upset(fromList(genedf_overlap[[3]][[2]]), 
            sets = names(genedf_overlap[[3]][[2]]), 
            order.by = "freq",
            mainbar.y.label = "# of overlapping genes", 
            sets.x.label = "# of genes",
            mb.ratio = c(0.5, 0.5))
p7 <- upset(fromList(genedf_overlap[[4]][[1]]), 
      sets = names(genedf_overlap[[4]][[1]]), 
      order.by = "freq",
      mainbar.y.label = "# of overlapping genes", 
      sets.x.label = "# of genes",
      mb.ratio = c(0.5, 0.5))
p8 <- upset(fromList(genedf_overlap[[4]][[2]]), 
            sets = names(genedf_overlap[[4]][[2]]), 
            order.by = "freq",
            mainbar.y.label = "# of overlapping genes", 
            sets.x.label = "# of genes",
            mb.ratio = c(0.5, 0.5))
p9 <- upset(fromList(genedf_overlap[[5]][[1]]), 
      sets = names(genedf_overlap[[5]][[1]]), 
      order.by = "freq",
      mainbar.y.label = "# of overlapping genes", 
      sets.x.label = "# of genes",
      mb.ratio = c(0.5, 0.5))
p10 <- upset(fromList(genedf_overlap[[5]][[2]]), 
            sets = names(genedf_overlap[[5]][[2]]), 
            order.by = "freq",
            mainbar.y.label = "# of overlapping genes", 
            sets.x.label = "# of genes",
            mb.ratio = c(0.5, 0.5))
p11 <- upset(fromList(genedf_overlap[[6]][[1]]), 
      sets = names(genedf_overlap[[6]][[1]]), 
      order.by = "freq",
      mainbar.y.label = "# of overlapping genes", 
      sets.x.label = "# of genes",
      mb.ratio = c(0.5, 0.5))
p12 <- upset(fromList(genedf_overlap[[6]][[2]]), 
            sets = names(genedf_overlap[[6]][[2]]), 
            order.by = "freq",
            mainbar.y.label = "# of overlapping genes", 
            sets.x.label = "# of genes",
            mb.ratio = c(0.5, 0.5))


pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/upset_geneoverlap_bulkMegasetandROSMAP.pdf"), 
    width = 5, height = 3.5)
p1 
grid.text("All AD vs Con, DFC, lower in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p2
grid.text("All AD vs Con, DFC, higher in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p3
grid.text("early AD vs Con, DFC, lower in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p4
grid.text("early AD vs Con, DFC, higher in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p5
grid.text("late AD vs early, DFC, lower in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p6
grid.text("late AD vs early, DFC, higher in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p7
grid.text("All AD vs Con, MTG, lower in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p8
grid.text("All AD vs Con, MTG, higher in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p9
grid.text("early AD vs Con, MTG, lower in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p10
grid.text("early AD vs Con, MTG, higher in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p11
grid.text("late AD vs early, MTG, lower in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
p12
grid.text("late AD vs early, MTG, higher in severe",x = 0.65, y=0.98, gp=gpar(fontsize=8))
dev.off()


#######
# Which modules are implicated in APOE comparison?
#####

apoe <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_APOE_DFC/shared_output_table.csv"))
#   mod Celltype_Gabitto_2024 
# 1  85                   Sst              
# 2 149           Endothelial              
# 3 744                   Vip   

# 85 and 149 are ribosomal modules
# Mod 744 is interesting
# > mods[[744]]
# [1] "RORA"   "SESN3"  "FUT9"   "LRRN1"  "MGAT4C"

###########
# Investigate modules shared between SCZ and AD
###########

title_vec <- c("Con vs all AD (DFC)",
             "Early AD vs Con (DFC)",
             "Late vs Early AD (DFC)",
             "Con vs all AD (MTG)",
             "Early AD vs Con (MTG)",
             "Late vs Early AD (MTG)",
             "Con vs all AD (ROSMAP, DFC)",
             "Early AD vs Con (ROSMAP, DFC)",
             "Late vs Early AD (ROSMAP, DFC)",
             "Con vs all AD (ROSMAP, MTG)",
             "Early AD vs Con (ROSMAP, MTG)",
             "Late vs Early AD (ROSMAP, MTG)",
             "APOE 4/4 vs 3/3 (DFC)",
             "SCZ vs con (DFC)")

dat_vec <- c("gabitto_vs_liu_DFC_AllADVsCon",    
             "gabitto_vs_liu_DFC_earlyVsCon",        
             "gabitto_vs_liu_DFC_lateVsEarly",        
             "gabitto_vs_liu_MTG_AllADVsCon",         
             "gabitto_vs_liu_MTG_earlyVsCon",   
             "gabitto_vs_liu_MTG_lateVsEarly",        
             "gabitto_vs_liu_DFC_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_DFC_earlyVsCon_ROSMAP", 
             "gabitto_vs_liu_DFC_lateVsEarly_ROSMAP",
             "gabitto_vs_liu_MTG_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_MTG_earlyVsCon_ROSMAP", 
             "gabitto_vs_liu_MTG_lateVsEarly_ROSMAP",
             "gabitto_vs_liu_APOE_DFC",
             "brainSCOPE_CMC_vs_SZBD") 

label_vec <- c("Shared in Gabitto and Liu DFC (con vs all AD)",    
               "Shared in Gabitto and Liu DFC (early AD vs con)",        
               "Shared in Gabitto and Liu DFC (late vs early AD)",        
               "Shared in Gabitto and Liu MTG (con vs all AD)",         
               "Shared in Gabitto and Liu MTG (early AD vs con)",   
               "Shared in Gabitto and Liu MTG (late vs early AD)",        
               "Shared in Gabitto and Liu DFC (con vs all AD, ROSMAP mods)",    
               "Shared in Gabitto and Liu DFC (early AD vs con, ROSMAP mods)", 
               "Shared in Gabitto and Liu DFC (late vs early AD, ROSMAP mods)",
               "Shared in Gabitto and Liu MTG (con vs all AD, ROSMAP mods)",    
               "Shared in Gabitto and Liu MTG (early AD vs con, ROSMAP mods)", 
               "Shared in Gabitto and Liu MTG (late vs early AD, ROSMAP mods)",
               "Shared in Gabitto and Liu DFC (APOE 4/4 vs 3/3)",
               "Shared in brainSCOPE CMC and SZBD DFC (SCZ vs con)") 

dat_list <- mapply(\(x, y){
  out <- fread(file.path(dir1, x, "shared_output_table.csv"), data.table = F)
  return(out)
}, dat_vec, title_vec, SIMPLIFY = F)
names(dat_list) <- label_vec

dat_list_neg <- lapply(dat_list, \(x) x[x$Direction == -1, ])
dat_list_pos <- lapply(dat_list, \(x) x[x$Direction == 1, ])

scz_mods_neg <- dat_list_neg[[14]]$mod # 23 mods
ad_mods_neg <- lapply(dat_list_neg[-14], \(x) x$mod) |> unlist() |> unique()
scz_mods_neg[scz_mods_neg %in% ad_mods_neg] # which scz mods (neg) are also implicated in AD?
#  [1]  34  78 100 105 132 187 274 298 316 324 329 329 382 405 489 489 489 489 507
# [20] 910
# Analyze # of overlapping mods per comparison:
lapply(dat_list_neg, \(x) sum(scz_mods_neg %in% unique(x$mod))) |> unlist()
            # Shared in Gabitto and Liu DFC (con vs all AD) 
            #                                                15 
            #   Shared in Gabitto and Liu DFC (early AD vs con) 
            #                                                 8 
            #  Shared in Gabitto and Liu DFC (late vs early AD) 
            #                                                 1 
            #     Shared in Gabitto and Liu MTG (con vs all AD) 
            #                                                17 
            #   Shared in Gabitto and Liu MTG (early AD vs con) 
            #                                                 6 
            #  Shared in Gabitto and Liu MTG (late vs early AD) 
            #                                                19 
#    Shared in Gabitto and Liu DFC (con vs all AD, ROSMAP mods) 
#                                                            14 
#  Shared in Gabitto and Liu DFC (early AD vs con, ROSMAP mods) 
#                                                             1 
# Shared in Gabitto and Liu DFC (late vs early AD, ROSMAP mods) 
#                                                             0 
#    Shared in Gabitto and Liu MTG (con vs all AD, ROSMAP mods) 
#                                                            11 
#  Shared in Gabitto and Liu MTG (early AD vs con, ROSMAP mods) 
#                                                             4 
# Shared in Gabitto and Liu MTG (late vs early AD, ROSMAP mods) 
#                                                             8 
#               Shared in Gabitto and Liu DFC (APOE 4/4 vs 3/3) 
#                                                             0 
#            Shared in brainSCOPE CMC and SZBD DFC (SCZ vs con) 
#                                                            23 



scz_mods_pos <- dat_list_pos[[14]]$mod # 3 mods
ad_mods_pos <- lapply(dat_list_pos[-14], \(x) x$mod) |> unlist() |> unique()
scz_mods_pos[scz_mods_pos %in% ad_mods_pos] # which scz mods (neg) are also implicated in AD?
# [1] 149
lapply(dat_list_pos, \(x) sum(scz_mods_pos %in% unique(x$mod))) |> unlist()
# Only 1 overlap (APOE 4/4 vs 3/3)
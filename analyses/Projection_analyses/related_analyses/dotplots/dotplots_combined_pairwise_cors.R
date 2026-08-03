library(tidyverse)
library(data.table)
library(qs)
library(ggpubr)
library(ComplexHeatmap)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/calculate_rand_euclidean_distances.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/dotplots/fxns.R"))

# Plot pairwise correlations between dotplot data.
# First, load all data:
dir1 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/")
dot_list <- list(
  file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_conVsAllAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_DFC.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_earlyADVscon.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_earlyADVscon_DFC.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_lateVsEarlyAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_lateVsEarlyAD_DFC.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_DFC_APOE", "dcopa_scorecard_summary_APOE44vs33.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_conVsAllAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_MTG.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_earlyADVscon.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_earlyADVscon_MTG.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_lateVsEarlyAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_lateVsEarlyAD_MTG.csv"),
  file.path(dir1, "brainSCOPE_CMC_SZBD_shared", "dcopa_scorecard_summary_SCZvsCon.csv")
) |>
  lapply(\(x) fread(x, data.table = F))
typevec_hi <- c("All AD vs Con DFC (n = 3)",
             "All AD vs Con ROSMAP DFC (n = 2)",
             "Early AD vs Con DFC (n = 2)",
             "Early AD vs Con ROSMAP DFC (n = 0)",
             "Late vs Early AD DFC (n = 0)",
             "Late vs Early AD ROSMAP, DFC (n = 4)",
             "APOE 4/4 vs 3/3 DFC (n = 2)",
             "All AD vs Con MTG (n = 3)",
             "All AD vs Con ROSMAP MTG (n = 0)",
             "Early AD vs Con MTG (n = 2)",
             "Early AD vs Con ROSMAP, MTG (n = 3)",
             "Late vs Early AD MTG (n = 2)",   
             "Late vs Early AD ROSMAP MTG (n = 2)",
             "SCZ vs Con DFC (n = 3)")

typevec_lo <- c("All AD vs Con DFC (n = 83)",
             "All AD vs Con ROSMAP DFC (n = 136)",
             "Early AD vs Con DFC (n = 15)",
             "Early AD vs Con ROSMAP DFC (n = 23)",
             "Late vs Early AD DFC (n = 8)",
             "Late vs Early AD ROSMAP, DFC (n = 16)",
             "APOE 4/4 vs 3/3 DFC (n = 1)",
             "All AD vs Con MTG (n = 135)",
             "All AD vs Con ROSMAP MTG (n = 186)",
             "Early AD vs Con MTG (n = 22)",
             "Early AD vs Con ROSMAP, MTG (n = 42)",
             "Late vs Early AD MTG (n = 119)",   
             "Late vs Early AD ROSMAP MTG (n = 139)",
             "SCZ vs Con DFC (n = 19)")

lapply(dot_list, \(x){
  x$type[grepl("Lower", x$type)] |> unique()
}) |> unlist()

dot_list <- mapply(\(l, name){
  l$comp <- name
  return(l)
}, dot_list, typevec, SIMPLIFY = F)

# Align SCZ celltypes to others
dot_list[[14]] <- dot_list[[14]] |> 
  mutate(Celltype = case_match(Celltype,
         "Astro" ~ "Astrocyte",
         "Oligo" ~ "Oligodendrocyte",
         "Micro" ~ "Microglia-PVM",
         "Endo" ~ "Endothelial",
         .default = Celltype)) |>
  filter(Celltype %in% dot_list[[1]]$Celltype)

# Ensure that celltypes are in same order
dot_list <- lapply(dot_list, \(x){
  x <- x |> mutate(Celltype = factor(Celltype,levels = unique(dot_list[[1]]$Celltype))) |>
    arrange(Celltype)
  return(x)
})

# Calculate correlation matrices for higher/lower in more severe
himat <- lapply(dot_list, \(x){
  x |> filter(grepl("Higher", type)) |> 
    mutate(num_sig = replace_na(num_sig, 0)) |>
    pull(num_sig)
}) |> do.call(what = "cbind") |>
  cor()
rownames(himat) <- typevec_hi
colnames(himat) <- typevec_hi

lomat <- lapply(dot_list, \(x){
  x |> filter(grepl("Lower", type)) |> 
    mutate(num_sig = replace_na(num_sig, 0)) |>
    pull(num_sig)
}) |> do.call(what = "cbind") |>
  cor()
rownames(lomat) <- typevec_lo
colnames(lomat) <- typevec_lo

himat[is.na(himat)] <- 0
lomat[is.na(lomat)] <- 0

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

phi <- Heatmap(himat,
               column_title = "Correlation between # of modules significantly\nupregulated in more severe cases",
               show_column_names = F,
               heatmap_legend_param = list(
                 direction = "horizontal",
                 title_gp = gpar(fontface = "plain"),
                 title_position = "topcenter"
               ),
               name = "Correlation",
               col = col_fun,    
               )
pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/dotplot_cor_higher.pdf"), width = 7.5, height = 4.5)
draw(phi, 
     heatmap_legend_side = "bottom",
     padding = unit(c(10, 30, 10, 30), "mm"))
dev.off()

plo <- Heatmap(lomat,
               column_title = "Correlation between # of modules significantly\ndownregulated in more severe cases",
               show_column_names = F,
               heatmap_legend_param = list(
                 direction = "horizontal",
                 title_gp = gpar(fontface = "plain"),
                 title_position = "topcenter"
               ),
               name = "Correlation",
               col = col_fun,    
               )
pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/summary_plots/dotplot_cor_lower.pdf"), width = 7.5, height = 4.5)
draw(plo, 
     heatmap_legend_side = "bottom",
     padding = unit(c(10, 30, 10, 30), "mm"))
dev.off()







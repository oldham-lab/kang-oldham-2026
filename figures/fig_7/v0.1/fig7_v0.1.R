library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)

version <- "v0.1"

##########
# Panel A-B
# Dotplots (Control vs all AD)
##########
# Find overlap between all AllADvsCon comparisons
file1 <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared"), full.names = T)
file1 <- file1[c(3:4, 10:11)]
flist <- lapply(file1, \(x) fread(data.table = F, file = file.path(x, "shared_output_table.csv")))
flistjoin <- Reduce(\(x, y){
  a1 <- x |> select(c("mod", "Celltype_Gabitto_2024", "Direction", "Consistency")) 
  a2 <- y |> select(c("mod", "Celltype_Gabitto_2024", "Direction", "Consistency"))
  out <- dplyr::inner_join(a1, a2, by = dplyr::join_by(mod, Celltype_Gabitto_2024, Direction, Consistency), relationship = "many-to-many") 
  return(out)
}, flist) |>
  filter(!duplicated(paste0(mod, Celltype_Gabitto_2024, Direction, Consistency))) |>
    mutate(Direction = factor(Direction, levels = c(1, -1))) |> 
    rename("Celltype" = "Celltype_Gabitto_2024") |>
    mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) |>
    select(mod, Celltype, Direction) |>
    group_by(Celltype, Direction, .drop = FALSE) |>
    summarise(num_sig = n(), .groups = "drop") |>
    left_join(class_info, by = join_by(Celltype == Subclass)) |>
    mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) |>
    arrange(Celltype) |>
    mutate(num_sig = case_match(num_sig, 
                                0 ~ NA,
                                .default = num_sig
                                ),
           comp = "Combined") 

## Plot dotplot
# For celltype ordering:
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

# dir1 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/")
# dot_list <- list(
#   file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_conVsAllAD.csv"),
#   file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_DFC.csv"),
#   file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_conVsAllAD.csv"),
#   file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_MTG.csv")
# ) |>
#   lapply(\(x){
#     fread(x, data.table = F) |> 
#       mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) 
#   })
# typevec_hi <- c("DFC (n = 3)",
#                 "ROSMAP DFC (n = 2)",
#                 "MTG (n = 3)",
#                 "ROSMAP MTG (n = 0)")

# typevec_lo <- c("DFC (n = 83)",
#                 "ROSMAP DFC (n = 136)",
#                 "MTG (n = 135)",
#                 "ROSMAP MTG (n = 186)"
#                 )

# df_hi <- lapply(seq_along(dot_list), \(x){
#   dot_list[[x]] |>
#     filter(grepl("Higher", type)) |>
#     mutate(comp = typevec_hi[x])
# }) |> do.call(what = "rbind")

# df_lo <- lapply(seq_along(dot_list), \(x){
#   dot_list[[x]] |>
#     filter(grepl("Lower", type)) |>
#     mutate(comp = typevec_lo[x])
# }) |> do.call(what = "rbind")

# cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

# phi <- df_hi |>
#     ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
#       theme_minimal() + 
#       geom_point(color = "black", pch = 21) +
#       theme(text = element_text(size = 6),
#             axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.2, 0, 0, 0, "cm")),
#             axis.text.y = element_text(size = 6), 
#             strip.text = element_text(size = 7),
#             legend.direction = "horizontal",    
#             legend.position = "bottom",
#             legend.box = "vertical",
#             legend.spacing.y = unit(4, "mm"),
#             legend.title = element_blank(),
#             legend.margin = margin(-0.5, 0, 0, 0, "cm"),
#             legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank()) +
#     #  guides(fill = guide_legend(ncol = 1),
#      #        size = guide_legend(ncol = 1)) +
#       labs(x = "", y = "") +
#       scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
#       scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) 
# ggsave(phi, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/panel_A1.pdf"), height = 2, width = 5)

# plo <- df_lo |>
#     ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
#       theme_minimal() + 
#       geom_point(color = "black", pch = 21) +
#       theme(text = element_text(size = 6),
#             axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.2, 0, 0, 0, "cm")),
#             axis.text.y = element_text(size = 6), 
#             strip.text = element_text(size = 7),
#             legend.direction = "horizontal",    
#             legend.position = "bottom",
#             legend.box = "vertical",
#             legend.spacing.y = unit(4, "mm"),
#             legend.title = element_blank(),
#             legend.margin = margin(-0.5, 0, 0, 0, "cm"),
#             legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank()) +
#     #  guides(fill = guide_legend(ncol = 1),
#      #        size = guide_legend(ncol = 1)) +
#       labs(x = "", y = "") +
#       scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
#       scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) 
# ggsave(plo, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/panel_A2.pdf"), height = 2.2, width = 5)

#############
# Try plotting single dataset dotplots
#############

# Load individual output tables

# con vs All DFC
olist <- list("g_dfc" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv")),
              "l_dfc" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_PFC/euclidean_distances/output_table_Subclass.csv")),
              "g_mtg" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv")),
              "l_mtg" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_MTC/euclidean_distances/output_table_Subclass.csv")),
              "g_r_dfc" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv")),
              "l_r_dfc" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv")),
              "g_r_mtg" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv")),
              "l_r_mtg" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv"))
              ) |> 
  lapply(\(x){
    x |> dplyr::filter(sig_FDR, Consistency %in% c(0, 1))
  })

# Reassign Liu celltype names to Gabitto
olist[c(2,4,6,8)] <- olist[c(2,4,6,8)] |>
  lapply(\(x){
    x |> 
      mutate(Celltype = case_match(Celltype,
        c("SMC", "VLMC", "End", "Per") ~ "Endothelial",
        c("Exc L4-5 IT-2", "Exc L3-4 IT","Exc L4-5 IT-1") ~ "L4 IT",
        "Exc L5 ET" ~ "L5 ET",
        c("Exc L4-5 IT-2", "Exc L4-5 IT-1","Exc L3-5 IT", "Exc L5-6 IT") ~ "L5 IT",
        "Exc L5/6 NP" ~ "L5/6 NP",
        "Inh LAMP5" ~ "Lamp5",
        "Inh PVALB" ~ "Pvalb",
        "Inh SST" ~ "Sst",
        "Exc L5-6 IT" ~ "L6 IT",
        "Exc L5/6 IT Car3" ~ "L6 IT Car3",
        "Exc L6 CT" ~ "L6 CT",
        "Inh PAX6" ~ "Pax6",
        "Ast" ~ "Astrocyte",
        "Inh VIP" ~ "Vip",
        "Exc L6b" ~ "L6b",
        "Exc L2-3 IT" ~ "L2/3 IT",
        "Exc L6 IT" ~ "L6 IT",
        "Mic" ~ "Microglia-PVM",
        "Oli" ~ "Oligodendrocyte",
        .default = Celltype
      )) |>
      filter(!Celltype %in% c("T", "Fib", "Exc EC"))
  })


# Create dotplots
olist_count <- lapply(olist, \(x){
  x |> 
    mutate(factor(Direction, levels = c(1, -1))) |>
    mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) |>
    select(mod, Celltype, Direction) |>
    group_by(Celltype, Direction, .drop = FALSE) |>
    summarise(num_sig = n(), .groups = "drop") |>
    left_join(class_info, by = join_by(Celltype == Subclass)) |>
    mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) |>
    arrange(Celltype) |>
    mutate(num_sig = case_match(num_sig, 
                                0 ~ NA,
                                .default = num_sig
                                )) 
})

typevec_hi <- c("Gabitto DFC",
                "Liu DFC",
                "Gabitto MTG",
                "Liu MTG",
                "Gabitto DFC (ROSMAP mods)",
                "Liu DFC (ROSMAP mods)",
                "Gabitto MTG (ROSMAP mods)",
                "Liu MTG (ROSMAP mods)")

typevec_lo <-  c("Gabitto DFC",
                "Liu DFC",
                "Gabitto MTG",
                "Liu MTG",
                "Gabitto DFC (ROSMAP mods)",
                "Liu DFC (ROSMAP mods)",
                "Gabitto MTG (ROSMAP mods)",
                "Liu MTG (ROSMAP mods)",
                "Combined")

ophi <- lapply(seq_along(olist_count), \(i){
  x <- olist_count[[i]]
  out <- x |> 
    filter(Direction == 1) |>
    mutate(comp = typevec_hi[i])
}) |> do.call(what = "rbind") |>
  mutate(comp = factor(comp, levels = rev(typevec_hi)))

oplo <- lapply(seq_along(olist_count), \(i){
  x <- olist_count[[i]]
  out <- x |> 
    filter(Direction == -1) |>
    mutate(comp = typevec_lo[i])
}) |> do.call(what = "rbind") 
oplo <- rbind(oplo, flistjoin |> filter(Direction == -1)) |>
  mutate(comp = factor(comp, levels = rev(typevec_lo)))

phi <- ophi |>
    ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 7), 
            legend.direction = "horizontal",    
            legend.position = "bottom",
            legend.box = "vertical",
            legend.spacing.y = unit(4, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(size = 7, margin = margin(0, 0, 0, -0.03, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
    #  guides(fill = guide_legend(ncol = 1),
     #        size = guide_legend(ncol = 1)) +
      labs(x = "", y = "") +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) 
ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_A.pdf"), height = 2.7, width = 5.5)

plo <- oplo |>
    ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 7), 
            legend.direction = "horizontal",    
            legend.position = "bottom",
            legend.box = "vertical",
            legend.spacing.y = unit(4, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(size = 7, margin = margin(0, 0, 0, -0.03, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
    #  guides(fill = guide_legend(ncol = 1),
     #        size = guide_legend(ncol = 1)) +
      labs(x = "", y = "") +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) 
ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7"), version, "panel_B.pdf"), height = 3, width = 5.5)

########
# Panel C
# Calculate pairwise hypergeometric p-value
########
nmod = 1023 * 24 * 2

# Fisher test (hypergeometric) function
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

# Find pairwise hypergeometric p-values
pmat <- matrix(nrow = length(olist), ncol = length(olist))
for(r in 1:nrow(pmat)){
  for(c in 1:ncol(pmat)){
    dat_out <- dplyr::inner_join(olist[[r]], olist[[c]], by = dplyr::join_by(mod, Celltype, Direction, Consistency), relationship = "many-to-many") |>
      dplyr::arrange(mod)
    pmat[r, c] <- fisherTest_modoverlap(unique(olist[[r]]$mod), # SEAAD2024
                                        unique(olist[[c]]$mod), # MIT
                                        length(unique(dat_out$mod)), # intersection between SEAAD2024 and MIT
                                        all = 1:nmod)
  }
}

# Log transform p-value mat
pmat_log <- -log10(pmat)
pmat_log[pmat_log == Inf] <- 308
rownames(pmat_log) <- typevec_hi
colnames(pmat_log) <- typevec_hi
sig_mat = matrix(ifelse(pmat_log > -log10(0.05/(nmod)), "*", ""), length(olist), length(olist))

p <- Heatmap(pmat_log,
             name = "-log(p-value)",
             cell_fun = function(j, i, x, y, w, h, col) {
               if(sig_mat[i, j] != "") {
                 grid.text(sig_mat[i, j], x, y)
               }
             }
)

pdf(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7"), version, "panel_C.pdf"), width = 6, height = 5)
draw(p)
dev.off()

##########
# Panel D - Upset plot
##########


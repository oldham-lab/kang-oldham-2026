library(tidyverse)
library(data.table)
library(qs)
library(ggpubr)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/calculate_rand_euclidean_distances.R"))

### Generate summary dotplot for sig mods that overlap between SEAAD2024 and MIT

# Function:

create_shared_dotplot_summary_AD <- function(slist,
                                             caption,
                                             splits,
                                             file_suffix,
                                             sn_anno_subclass,
                                             ct_order,
                                             rev = F,
                                             save_dir1){
  

  #sigmod_count <- lapply(slist, \(x) length(unique(x$mod))) # Count of significant mods
  sigmod_count <- c(slist |> filter(type=="pos") |> pull(mod) |> unique() |> length(),
                    slist |> filter(type=="neg") |> pull(mod) |> unique() |> length())

  # hi_sum <- slist[[1]] |> 
  #   group_by(subclass) |> 
  #   summarise("num_sig" = n()) |>
  #   right_join(sn_anno_subclass, by = join_by(subclass == Subclass)) |>
  #   mutate(num_sig = replace_na(num_sig, 0),
  #          type =  paste0("Higher in ", splits[1], " vs ", splits[2], "\n(n = ",sigmod_count[[1]], ")"))

  # lo_sum <- slist[[2]] |> 
  #   group_by(subclass) |> 
  #   summarise("num_sig" = n()) |>
  #   right_join(sn_anno_subclass, by = join_by(subclass == Subclass)) |>
  #   mutate(num_sig = replace_na(num_sig, 0),
  #          type =  paste0("Lower in ", splits[1], " vs ", splits[2],  "\n(n = ",sigmod_count[[2]], ")"))

  hi_sum <- slist |> 
    filter(type == "pos") |>
    group_by(subclass) |> 
    summarise("num_sig" = n()) |>
    right_join(sn_anno_subclass, by = join_by(subclass == Subclass)) |>
    mutate(num_sig = replace_na(num_sig, 0),
          type =  paste0("Higher in ", splits[1], " vs ", splits[2], "\n(n = ",sigmod_count[[1]], ")")) |>
    arrange(factor(Subclass_fixed, levels = ct_order))

  lo_sum <- slist |> 
    filter(type == "neg") |>
    group_by(subclass) |> 
    summarise("num_sig" = n()) |>
    right_join(sn_anno_subclass, by = join_by(subclass == Subclass)) |>
    mutate(num_sig = replace_na(num_sig, 0),
          type =  paste0("Lower in ", splits[1], " vs ", splits[2],  "\n(n = ",sigmod_count[[2]], ")")) |>
    arrange(factor(Subclass_fixed, levels = ct_order))

  cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

  if(rev){
    allsum <- rbind(lo_sum, hi_sum)
  } else {
    allsum <- rbind(hi_sum, lo_sum)
  }


  psumall <- allsum |>
    mutate(pcut = "",
          type = factor(type, levels = rev(unique(type))),
          subclass = factor(subclass, levels = unique(subclass))) |>
    ggplot(aes(x = subclass, y = type, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.2, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 6), 
            strip.text = element_text(size = 7),
            legend.position = "bottom",
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(x = "", y = "",
          title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
          subtitle = "1023 total modules, FDR cutoff",
          caption = caption) +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) 

  ggsave(psumall, file = paste0(save_dir1, "/dcopa_scorecard_summary_", file_suffix, ".pdf"), width = 5.5, height = 2.2, bg = "white", limitsize=F)
  ggsave(psumall, file = paste0(save_dir1, "/dcopa_scorecard_summary_", file_suffix, ".png"), width = 5.5, height = 2.2, bg = "white", limitsize=F)
}

sn_anno_subclass <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |> 
  arrange(Subclass) |>
  mutate(Subclass_fixed = c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip")) 

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC")

# Load sig modules (conVAll)
save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared/")
if(!dir.exists(save_dir1)){dir.create(save_dir1)}

overlap1 <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/seaad2024_mit_dcopa_overlap_PFC.qs"))
overlap1 <- lapply(overlap1, \(x){
  out <- x |> 
    rename("type" = "type.x", "subclass" = "subclass.x")
}) # use Gabitto subclass names


create_shared_dotplot_summary_AD(slist = overlap1[[1]],
                                 caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); control vs all AD samples",
                                 splits = c("all AD", "con"),
                                 file_suffix = "conVsAllAD",
                                 sn_anno_subclass = sn_anno_subclass,
                                 ct_order = ct_order,
                                 rev = F,
                                 save_dir1 = save_dir1)

create_shared_dotplot_summary_AD(slist = overlap1[[2]],
                                 caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); control vs early AD samples",
                                 splits = c("con", "early AD"),
                                 file_suffix = "conVsEarly",
                                 sn_anno_subclass = sn_anno_subclass,
                                 ct_order = ct_order,
                                 rev = T,
                                 save_dir1 = save_dir1)

create_shared_dotplot_summary_AD(slist = overlap1[[3]],
                                 caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); early vs late AD samples",
                                 splits = c("late AD", "early AD"),
                                 file_suffix = "earlyVsLate",
                                 sn_anno_subclass = sn_anno_subclass,
                                 ct_order = ct_order,
                                 rev = F,
                                 save_dir1 = save_dir1)

# How many modules are shared in common across shared AD and shared SCZ modules?
overlap_scz <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/brainscope_cmc_szbd_dcopa_overlap.qs"))

length(unique(overlap1[[1]]$mod))
# 89 unique mods for AD conVAll (86 lower in AD, 3 higher in AD)
length(unique(overlap_scz[[1]]$mod))
# 21 unique mods for SCZ (18 lower in SCZ, 3 higher in SCZ)
sum(unique(overlap_scz[[1]]$mod) %in% unique(overlap1[[1]]$mod))
# 13 mods shared between the two (all lower in AD)

all_subclasses <- unique(c(overlap1[[1]]$subclass, overlap_scz[[1]]$subclass))
overlap_join <- full_join(overlap1[[1]][,1:3], overlap_scz[[1]][,1:3], by = join_by(mod))
overlap_join1 <- overlap_join |> filter(!is.na(subclass) & !is.na(subclass.x)) |>
  filter(subclass == subclass.x)
# 7 unique modules that move in same direction and same subclass:
# > overlap_join1
#    mod subclass type subclass.x type.x
# 1   34    L4 IT  neg      L4 IT    neg
# 2  100    L4 IT  neg      L4 IT    neg
# 3  100    L4 IT  neg      L4 IT    neg
# 4  132    L4 IT  neg      L4 IT    neg
# 5  132    L4 IT  neg      L4 IT    neg
# 6  324    L4 IT  neg      L4 IT    neg
# 7  324    L4 IT  neg      L4 IT    neg
# 8  329    L4 IT  neg      L4 IT    neg
# 9  329    L4 IT  neg      L4 IT    neg
# 10 405    L4 IT  neg      L4 IT    neg
# 11 489    L5 IT  neg      L5 IT    neg
# 12 489    L5 IT  neg      L5 IT    neg
# 13 489  L5/6 NP  neg    L5/6 NP    neg

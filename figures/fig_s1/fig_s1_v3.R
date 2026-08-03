# v3: add panels for DFC

library(qs)
library(data.table)
library(tidyverse)
library(ggrepel)
library(cowplot)

# Scatterplots comparing the proportions of DFC/MTG nuclei assigned to each shared 
# cell class (a), subclass (b), or supertype (c) in Jorstad et al. and Gabitto et al.

#####
# DFC
#####
jor <- fread("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_DFC.csv", data.table = FALSE) |>
  mutate(Class = case_when(
    Cell_Type %in% c("Astro", "Endo", "Micro/PVM", "Oligo", "OPC", "VLMC") ~ "Non-neuronal",
    Cell_Type %in% c("Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip") ~ "GABAergic",
    TRUE ~ "Glutamatergic"
  )) |>
  rename("Subclass" = "Cell_Type",
         "Supertype" = "Cluster")
gab <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  mutate(
    Class = case_match(Class,
      "Neuronal: GABAergic" ~ "GABAergic",
      "Neuronal: Glutamatergic" ~ "Glutamatergic",
      "Non-neuronal and Non-neural" ~ "Non-neuronal"
    ),
    Subclass = case_match(Subclass,
      "Astrocyte" ~ "Astro",
      "Endothelial" ~ "Endo",
      "Microglia-PVM" ~ "Micro/PVM",
      "Oligodendrocyte" ~ "Oligo",
      .default = Subclass),
    Supertype = case_match(Supertype,
      "Lamp5_Lhx6_1" ~ "Lamp5 Lhx6_1",
      .default = Supertype),
    Supertype = gsub("Micro-PVM", "Micro/PVM", Supertype)
  )

# Calculate proportions
plotdflist <- lapply(c("Class", "Subclass", "Supertype"), \(x){
  jordf <- jor |> group_by_at(x) |>
    summarise(count = n()) |>
    mutate(prop_jor = count/sum(count))
  gabdf <- gab |> group_by_at(x) |>
    summarise(count = n()) |>
    mutate(prop_gab = count/sum(count))
  outdf <- full_join(jordf[, c(1,3)], gabdf[, c(1,3)]) |>
    mutate(
      prop_jor = case_when(
        is.na(prop_jor) ~ 0,
        .default = prop_jor
      ),
      prop_gab = case_when(
        is.na(prop_gab) ~ 0,
        .default = prop_gab
      ),
      type = x)
  colnames(outdf)[1] <- "Label"
  return(as.data.frame(outdf))
})

r_vec <- data.frame("type" = c("Class", "Subclass", "Supertype"),
                    "lab" = lapply(plotdflist, \(x) paste0("r = ", sprintf("%.2g", cor(x$prop_jor, x$prop_gab)))) |> unlist(),
                    "x_pos" = c(0.5, 0.2, 0.2))

# Plot figure
pall <- ggplot(plotdflist |> do.call(what = "rbind"), aes(x = prop_jor, y = prop_gab, label = Label)) +
  theme_classic() +
  geom_point() +
  geom_abline(linetype = "dashed", alpha = 0.5) +
  geom_label_repel(fill = "white",
                   min.segment.length = 0) +
  geom_text(data = r_vec, aes(label = lab, x = x_pos, y = 0)) +
  labs(x = "Proportion of nuclei\n(Jorstad et al., 2023, DFC)",
       y = "Proportion of nuclei\n(Gabitto et al., 2024, DFC)") +
  theme(text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        strip.text = element_text(size = 12)) +
  facet_wrap(~type, nrow = 1, scales = "free")
ggsave(pall, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s1/v3/panel_AtoC.pdf")), height = 3.5, width = 10)
ggsave(pall, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s1/v3/panel_AtoC.png")), height = 3.5, width = 10)


#####
# MTG
#####
jor <- fread("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_MTG.csv", data.table = FALSE) |>
  mutate(Class = case_when(
    Cell_Type %in% c("Astro", "Endo", "Micro/PVM", "Oligo", "OPC", "VLMC") ~ "Non-neuronal",
    Cell_Type %in% c("Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip") ~ "GABAergic",
    TRUE ~ "Glutamatergic"
  )) |>
  rename("Subclass" = "Cell_Type",
         "Supertype" = "Cluster")
gab <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  mutate(
    Class = case_match(Class,
      "Neuronal: GABAergic" ~ "GABAergic",
      "Neuronal: Glutamatergic" ~ "Glutamatergic",
      "Non-neuronal and Non-neural" ~ "Non-neuronal"
    ),
    Subclass = case_match(Subclass,
      "Astrocyte" ~ "Astro",
      "Endothelial" ~ "Endo",
      "Microglia-PVM" ~ "Micro/PVM",
      "Oligodendrocyte" ~ "Oligo",
      .default = Subclass),
    Supertype = case_match(Supertype,
      "Lamp5_Lhx6_1" ~ "Lamp5 Lhx6_1",
      .default = Supertype),
    Supertype = gsub("Micro-PVM", "Micro/PVM", Supertype)
  )

# Calculate proportions
plotdflist <- lapply(c("Class", "Subclass", "Supertype"), \(x){
  jordf <- jor |> group_by_at(x) |>
    summarise(count = n()) |>
    mutate(prop_jor = count/sum(count))
  gabdf <- gab |> group_by_at(x) |>
    summarise(count = n()) |>
    mutate(prop_gab = count/sum(count))
  outdf <- full_join(jordf[, c(1,3)], gabdf[, c(1,3)]) |>
    mutate(
      prop_jor = case_when(
        is.na(prop_jor) ~ 0,
        .default = prop_jor
      ),
      prop_gab = case_when(
        is.na(prop_gab) ~ 0,
        .default = prop_gab
      ),
      type = x)
  colnames(outdf)[1] <- "Label"
  return(as.data.frame(outdf))
})

r_vec <- data.frame("type" = c("Class", "Subclass", "Supertype"),
                    "lab" = lapply(plotdflist, \(x) paste0("r = ", sprintf("%.2g", cor(x$prop_jor, x$prop_gab)))) |> unlist(),
                    "x_pos" = c(0.5, 0.2, 0.2))

# Plot figure
pall <- ggplot(plotdflist |> do.call(what = "rbind"), aes(x = prop_jor, y = prop_gab, label = Label)) +
  theme_classic() +
  geom_point() +
  geom_abline(linetype = "dashed", alpha = 0.5) +
  geom_label_repel(fill = "white",
                   min.segment.length = 0) +
  geom_text(data = r_vec, aes(label = lab, x = x_pos, y = 0)) +
  labs(x = "Proportion of nuclei\n(Jorstad et al., 2023, MTG)",
       y = "Proportion of nuclei\n(Gabitto et al., 2024, MTG)") +
  theme(text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        strip.text = element_text(size = 12)) +
  facet_wrap(~type, nrow = 1, scales = "free")
ggsave(pall, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s1/v3/panel_DtoF.pdf")), height = 3.5, width = 10)
ggsave(pall, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s1/v3/panel_DtoF.png")), height = 3.5, width = 10)

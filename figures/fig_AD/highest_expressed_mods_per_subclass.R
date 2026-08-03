library(DESeq2)
library(qs)
library(data.table)
library(AnnotationHub)
library(tidyverse)
library(cowplot)
library(data.table)
options(bitmapType = 'cairo')

# Load SEA projections
con_ind <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/sn_proj_indices/log_native/indices_over_all_datasets_Subclass_Control.csv"))
ad_ind <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/sn_proj_indices/log_native/indices_over_all_datasets_Subclass_Alzheimers.csv"))

# For each subclass find # of mods where subclass is highest expressed
conmax <- apply(con_ind[ ,-25], 1, \(x) colnames(con_ind)[which.max(x)]) |> table() |> c() |> 
  stack() |>
  dplyr::mutate(type = "Con")
admax <- apply(ad_ind[ ,-25], 1, \(x) colnames(con_ind)[which.max(x)]) |> table() |> c() |> 
  stack() |>
  dplyr::mutate(type = "AD")
plotdf <- rbind(conmax, admax)
plotord <- plotdf |> 
  group_by(ind) |>
  summarise(mean_count = mean(values)) |>
  arrange(desc(mean_count))
plotdf$ind <- factor(plotdf$ind, levels = plotord$ind)

p <- ggplot(plotdf, aes(x = ind, y = values, fill = type)) +
  theme_classic() +
  geom_bar(position = "dodge", stat = "identity") +
  labs(x = "", y = "# of modules",
       title = "# of modules with highest expression in each subclass",
       subtitle = "Log native expression of bulk megaset modules in SEAAD2024") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.5),
        legend.title = element_blank(),
        legend.box.background = element_rect(color = "black", linewidth = 1),
        plot.title = element_text(face = "bold"))
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/count_per_subclass_of_mods_with_max_expr.pdf"), height = 3.5, width = 6)
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/count_per_subclass_of_mods_with_max_expr.png"), height = 3.5, width = 6)

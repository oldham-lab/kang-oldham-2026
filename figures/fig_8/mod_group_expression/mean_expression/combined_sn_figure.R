# Combined single-nucleus FEATURED figure: a summary table of modules + genes per
# module group, stacked above the per-shared-module single-nucleus boxplots
# (shared_module_featured_sn) and their companion Cliff's-delta stats heatmap.
# Sourcing shared_module_featured_figure.R rebuilds p_sn (boxplots) and p_sn_stats
# (heatmap) exactly as the standalone figures (it also re-saves them; harmless),
# and brings gene_tbl / group_overall into scope for the table.

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression/shared_module_featured_figure.R"))

library(patchwork)
library(ggpubr)

# Table: unique modules and genes per module group (group_overall, collapsed over
# cell types) -- describes how the gene/module universe splits by AD/SCZ overlap.
tbl <- gene_tbl %>%
  mutate(Group = group_overall(mod)) %>%
  group_by(Group) %>%
  summarise(Modules = n_distinct(mod), Genes = n(), .groups = "drop") %>%
  arrange(Group)
tbl_disp <- bind_rows(
  mutate(tbl, Group = as.character(Group)),
  tibble(Group = "Total", Modules = sum(tbl$Modules), Genes = sum(tbl$Genes))
) %>%
  mutate(Modules = formatC(Modules, big.mark = ",", format = "d"),
         Genes   = formatC(Genes,   big.mark = ",", format = "d"))

p_tbl <- ggtexttable(tbl_disp, rows = NULL,
                     theme = ttheme("light", base_size = 18)) +
  labs(title = "How many modules and genes fall in each AD/SCZ-overlap group",
       caption = "Shared = DE in both AD & SCZ (same cell type & direction); Unshared = DE in one disease only; Neither = not DE.  Universe: 1,016 modules / 17,972 genes.") +
  theme(plot.title   = element_text(size = 17, face = "bold", hjust = 0.5),
        plot.caption = element_text(size = 11, face = "italic", hjust = 0.5),
        plot.margin  = margin(10, 10, 10, 10))

# Clearer, plain-language titles on the reused featured panels.
p_box <- p_sn + labs(
  title = "Shared-module genes vs other genes: expression by cell type (single-nucleus)",
  subtitle = "One panel = a module disrupted in both AD & SCZ; 'Module' = its genes, vs Unshared & Neither genes.") +
  theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 13))

p_stat <- p_sn_stats + labs(
  title = "Effect size (Cliff's delta): shared-module genes vs other genes (single-nucleus)") +
  theme(plot.title = element_text(size = 15))

# Stack: table (small), boxplots (large), heatmap (medium).
combined <- p_tbl / p_box / p_stat +
  plot_layout(heights = c(0.12, 1, 0.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

out <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression/mod_group_expression_combined_sn.pdf")
ggsave(out, combined, width = 13, height = 33, limitsize = FALSE)
cat("saved", out, "\n")

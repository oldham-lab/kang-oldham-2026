# Combined PC1-variance FEATURED figure (equivalent of the mean-expression
# combined single-nucleus figure): a summary table of modules + genes per module
# group, stacked above the per-shared-module featured figure (pc1_variance_featured,
# bulk). PC1 variance is bulk-only and its stats are shown in-plot as Wilcoxon
# brackets, so there is no separate companion stats panel. Sourcing
# pc1_variance_featured_figure.R rebuilds `p` exactly as the standalone figure (and
# re-saves it; harmless) and brings gene_tbl / group_overall into scope.

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance/pc1_variance_featured_figure.R"))

library(patchwork)
library(ggpubr)

# Table: unique modules and genes per module group (group_overall, collapsed over
# cell types).
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

# Clearer, plain-language title on the reused featured panel.
p_feat <- p + labs(
  title = "Shared-module genes vs other genes: how tightly co-expressed they are (bulk)",
  subtitle = "Per-gene PC1 variance explained (higher = more tightly co-expressed); one panel = a module disrupted in both AD & SCZ.") +
  theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 13))

combined <- p_tbl / p_feat +
  plot_layout(heights = c(0.10, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

out <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance/pc1_variance_combined_featured.pdf")
ggsave(out, combined, width = 13, height = 24, limitsize = FALSE)
cat("saved", out, "\n")

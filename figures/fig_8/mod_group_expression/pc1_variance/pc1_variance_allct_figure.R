# Figure + Wilcoxon: per-module PC1 % variance explained (bulk), compared across
# the Shared / Unshared / Neither module groups, COLLAPSED over cell type. The
# group is cell-type-AGNOSTIC (group_overall in ../mod_group_defs.R): a module is
# Shared if it is a shared (mod, ct, dir) triplet in ANY cell type, Unshared if DE
# in any cell type in only one disease, else Neither. Unlike pc1_variance_figure.R
# (which replicates each module across all 23 cell types and facets), here each
# module contributes ONE value to ONE group, shown in a single pooled panel.

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance/pc1_variance_defs.R"))

fig_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance")

# Per-module PVE (from pc1_variance_defs.R), one cell-type-agnostic group each.
all_long <- pve_tbl %>%
  filter(!is.na(pve)) %>%
  mutate(group = group_overall(mod), series = "Bulk", stratum = "All cell types") %>%
  filter(!is.na(group))
fwrite(all_long, file.path(fig_dir, "pc1_variance_allct_long.csv"))
cat("modules per group:\n"); print(table(all_long$group))

# ---------------------------------------------------------------------------
# Wilcoxon between groups (one pooled stratum) with Cliff's delta.
# ---------------------------------------------------------------------------
pairs <- list(c("Shared", "Neither"), c("Unshared", "Neither"), c("Shared", "Unshared"))

wilcox_one <- function(d, g1, g2) {
  x <- d$pve[d$group == g1]; y <- d$pve[d$group == g2]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  w <- wilcox.test(x, y); U <- unname(w$statistic)
  tibble(group1 = g1, group2 = g2, n1 = length(x), n2 = length(y),
         median1 = median(x), median2 = median(y), median_diff = median(x) - median(y),
         W = U, cliffs_delta = 2 * U / (length(x) * length(y)) - 1, p = w$p.value)
}

wilcox_res <- all_long %>%
  mutate(dataset = "Bulk", condition = "Bulk") %>%
  group_by(dataset, condition, stratum) %>%
  group_modify(~ map_dfr(pairs, function(pr) wilcox_one(.x, pr[1], pr[2]))) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH"))
fwrite(wilcox_res, file.path(fig_dir, "pc1_variance_allct_wilcoxon.csv"))
cat("\nWilcoxon tests:", nrow(wilcox_res), "| BH<0.05:", sum(wilcox_res$p_adj_BH < 0.05), "\n")

# ---------------------------------------------------------------------------
# Figure: single pooled panel, one box per group, Wilcoxon brackets above whiskers.
# ---------------------------------------------------------------------------
lay <- bracket_layout(signif_brackets(wilcox_res, levels(all_long$group)),
                      all_long, "stratum", "pve")

ncount <- all_long %>% group_by(stratum, group) %>%
  summarise(n = n_distinct(mod), .groups = "drop")
rng <- lay$ylim[2] - lay$ylim[1]
ncount$y <- lay$ylim[1] - 0.05 * rng
lay$ylim[1] <- lay$ylim[1] - 0.11 * rng

p <- ggplot(all_long, aes(x = group, y = pve, fill = series)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.3, width = 0.6) +
  geom_signif(data = lay$sig, manual = TRUE, inherit.aes = FALSE,
              textsize = 3.4, tip_length = 0.01, vjust = -0.3,
              aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y_position)) +
  geom_text(data = ncount, inherit.aes = FALSE,
            aes(x = group, y = y, label = paste0("n=", n)),
            size = 3.2, vjust = 1, color = "grey25") +
  coord_cartesian(ylim = lay$ylim) +
  facet_wrap(~ stratum) +
  labs(x = NULL, y = "PC1 variance explained (%, mean per-gene R²)",
       fill = "Dataset",
       title = "Per-module PC1 variance explained by module group (bulk)",
       subtitle = "Module = one point; cell-type-agnostic groups, pooled over all cell types",
       caption = "Brackets: Wilcoxon BH-adjusted p (*** <0.001, ** <0.01, * <0.05).") +
  theme_bw(base_size = 17) +
  theme(axis.text.x  = element_text(size = 15),
        axis.text.y  = element_text(size = 14),
        axis.title   = element_text(size = 16),
        strip.text   = element_text(size = 16),
        plot.title   = element_text(size = 19, face = "bold"),
        plot.subtitle = element_text(size = 16),
        legend.position = "none")

ggsave(file.path(fig_dir, "pc1_variance_allct.pdf"), p, width = 10, height = 8)
ggsave(file.path(fig_dir, "pc1_variance_allct.png"), p, width = 10, height = 8, dpi = 200)
cat("done\n")

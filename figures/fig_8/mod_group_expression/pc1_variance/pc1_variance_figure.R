# Figure + Wilcoxon: per-module PC1 % variance explained (bulk), compared across
# the Shared / Unshared / Neither module groups -- same cell-type-resolved
# grouping as the mean_expression analysis (see ../mod_group_defs.R).
#
# Metric (one value per module, bulk only): take the bulk expression of the
# module's genes (genes x samples), run PCA over samples (correlation PCA), and
# average, across the module's genes, the % of each gene's variance explained by
# PC1 (= cor(gene, PC1)^2). With correlation PCA this mean equals PC1's overall
# proportion of variance explained. Bulk is celltype-agnostic, so each module's
# value is replicated across all cell types and re-grouped by that cell type's
# membership, then faceted by cell type (mirrors the bulk mean-expression figure).

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance/pc1_variance_defs.R"))

fig_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance")

# Per-module PVE + per-gene R² come from pc1_variance_defs.R
fwrite(pve_tbl, file.path(fig_dir, "pc1_variance_per_module.csv"))
fwrite(r2_tbl,  file.path(fig_dir, "pc1_variance_per_gene.csv"))
cat("modules:", nrow(pve_tbl), "| with PVE:", sum(!is.na(pve_tbl$pve)),
    "| median PVE:", round(median(pve_tbl$pve, na.rm = TRUE), 1), "%\n")

# ---------------------------------------------------------------------------
# Replicate across cell types and assign the cell-type-resolved group
# ---------------------------------------------------------------------------
all_long <- pve_tbl %>%
  filter(!is.na(pve)) %>%
  crossing(celltype = expr_celltypes) %>%
  mutate(group = group_in_celltype(mod, celltype), series = "Bulk") %>%
  filter(!is.na(group))
fwrite(all_long, file.path(fig_dir, "pc1_variance_long.csv"))
cat("rows per group:\n"); print(table(all_long$group))

# ---------------------------------------------------------------------------
# Wilcoxon between groups (per cell type) with Cliff's delta. Run before the
# figure so it can annotate the brackets.
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
  group_by(dataset, condition, celltype) %>%
  group_modify(~ map_dfr(pairs, function(pr) wilcox_one(.x, pr[1], pr[2]))) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH")) %>%
  arrange(celltype, group1, group2)
fwrite(wilcox_res, file.path(fig_dir, "pc1_variance_wilcoxon.csv"))
cat("\nWilcoxon tests:", nrow(wilcox_res), "| BH<0.05:", sum(wilcox_res$p_adj_BH < 0.05), "\n")
cat("Shared>Neither (positive Cliff's delta):",
    with(subset(wilcox_res, group1 == "Shared" & group2 == "Neither"), sum(cliffs_delta > 0)),
    "/", nrow(subset(wilcox_res, group1 == "Shared" & group2 == "Neither")), "\n")

# ---------------------------------------------------------------------------
# Figure: boxplots faceted by cell type, with Wilcoxon brackets (bulk = one box
# per group, so brackets are unambiguous; mirrors the bulk mean-expression fig).
# ---------------------------------------------------------------------------
lay <- bracket_layout(signif_brackets(wilcox_res, levels(all_long$group)),
                      all_long, "celltype", "pve")   # brackets above whiskers, stacked

# Modules analyzed per group (unit = module here), labelled below the boxes.
ncount <- all_long %>% group_by(celltype, group) %>%
  summarise(n = n_distinct(mod), .groups = "drop")
rng <- lay$ylim[2] - lay$ylim[1]
ncount$y <- lay$ylim[1] - 0.05 * rng
lay$ylim[1] <- lay$ylim[1] - 0.11 * rng

p <- ggplot(all_long, aes(x = group, y = pve, fill = series)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.3) +
  geom_signif(data = lay$sig, manual = TRUE, inherit.aes = FALSE,
              textsize = 3.0, tip_length = 0.01, vjust = -0.3,
              aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y_position)) +
  geom_text(data = ncount, inherit.aes = FALSE,
            aes(x = group, y = y, label = paste0("n=", n)),
            size = 2.7, vjust = 1, color = "grey25") +
  coord_cartesian(ylim = lay$ylim) +
  facet_wrap(~ celltype) +
  labs(x = NULL, y = "PC1 variance explained (%, mean per-gene R²)",
       fill = "Dataset",
       title = "Per-module PC1 variance explained by module group (bulk)",
       subtitle = "Module = one point; cell-type-resolved groups, all cell types",
       caption = "Brackets: Wilcoxon BH-adjusted p (*** <0.001, ** <0.01, * <0.05).") +
  theme_bw(base_size = 17) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title   = element_text(size = 16),
        strip.text   = element_text(size = 15),
        plot.title   = element_text(size = 19, face = "bold"),
        plot.subtitle = element_text(size = 16),
        legend.title = element_text(size = 15),
        legend.text  = element_text(size = 14),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "pc1_variance_by_celltype.pdf"), p, width = 16, height = 16)
ggsave(file.path(fig_dir, "pc1_variance_by_celltype.png"), p, width = 16, height = 16, dpi = 200)
cat("done\n")

# Figure + Wilcoxon: mean gene expression across the three module groups,
# COLLAPSED over cell type (no per-cell-type facet). The group is cell-type-
# AGNOSTIC (group_overall in ../mod_group_defs.R): a module is Shared if it is a
# shared (mod, ct, dir) triplet in ANY cell type, Unshared if DE in any cell type
# in only one disease, else Neither. Single-nucleus expression is POOLED over all
# 23 cell types (every gene x cell type value is a point); bulk is one value per
# gene. Bulk and single-nucleus rendered as separate figures (each its own
# y-axis), mirroring mod_group_expression_figure.R.

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression/mean_expression_defs.R"))

fig_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression")

# Pooled long table: sc keeps every (gene, cell type) value; bulk is one per gene
# (NOT replicated across cell types). Group is cell-type-agnostic.
pooled <- bind_rows(
  sc_z %>%
    inner_join(gene_tbl, by = "gene") %>%
    mutate(group = group_overall(mod), series = paste(dataset, condition)),
  bulk_z %>%
    inner_join(gene_tbl, by = "gene") %>%
    mutate(dataset = "Bulk", condition = "Bulk", series = "Bulk",
           group = group_overall(mod))
) %>%
  filter(!is.na(group)) %>%
  mutate(stratum = "All cell types")

fwrite(filter(pooled, dataset != "Bulk"), file.path(fig_dir, "mod_group_expression_allct_long_sn.csv"))
fwrite(filter(pooled, dataset == "Bulk"), file.path(fig_dir, "mod_group_expression_allct_long_bulk.csv"))
cat("rows per group:\n"); print(table(pooled$group))

# ---------------------------------------------------------------------------
# Wilcoxon between groups per series (= dataset x condition), with Cliff's delta.
# (sc points are pooled over cell types, so n is large and p tiny -- Cliff's delta
# is the meaningful effect size; same convention as the per-cell-type figures.)
# ---------------------------------------------------------------------------
pairs <- list(c("Shared", "Neither"), c("Unshared", "Neither"), c("Shared", "Unshared"))

wilcox_one <- function(d, g1, g2) {
  x <- d$expr_scale[d$group == g1]; y <- d$expr_scale[d$group == g2]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  w <- wilcox.test(x, y); U <- unname(w$statistic)
  np <- as.numeric(length(x)) * as.numeric(length(y))   # numeric: avoid integer overflow at large n
  tibble(group1 = g1, group2 = g2, n1 = length(x), n2 = length(y),
         median1 = median(x), median2 = median(y), median_diff = median(x) - median(y),
         W = U, cliffs_delta = 2 * U / np - 1, p = w$p.value)
}

wilcox_res <- pooled %>%
  group_by(dataset, condition) %>%
  group_modify(~ map_dfr(pairs, function(pr) wilcox_one(.x, pr[1], pr[2]))) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH"),
         stratum = "All cell types") %>%
  arrange(dataset, condition, group1, group2)
fwrite(filter(wilcox_res, dataset != "Bulk"), file.path(fig_dir, "mod_group_expression_allct_wilcoxon_sn.csv"))
fwrite(filter(wilcox_res, dataset == "Bulk"), file.path(fig_dir, "mod_group_expression_allct_wilcoxon_bulk.csv"))
cat("\nWilcoxon tests:", nrow(wilcox_res), "| BH<0.05:", sum(wilcox_res$p_adj_BH < 0.05), "\n")

# ---------------------------------------------------------------------------
# Bulk figure: single pooled panel, one box per group, Wilcoxon brackets.
# ---------------------------------------------------------------------------
bulk <- filter(pooled, dataset == "Bulk")
cap  <- quantile(bulk$z, c(0.01, 0.97), na.rm = TRUE)        # winsorize for display
bulk <- mutate(bulk, z = pmin(pmax(z, cap[1]), cap[2]))
lay  <- bracket_layout(signif_brackets(filter(wilcox_res, dataset == "Bulk"),
                                       levels(pooled$group)), bulk, "stratum", "z")
ncount_b <- bulk %>% group_by(stratum, group) %>% summarise(n = n_distinct(gene), .groups = "drop")
rng <- lay$ylim[2] - lay$ylim[1]
ncount_b$y <- lay$ylim[1] - 0.05 * rng
lay$ylim[1] <- lay$ylim[1] - 0.11 * rng

p_bulk <- ggplot(bulk, aes(x = group, y = z, fill = series)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.3, width = 0.6) +
  geom_signif(data = lay$sig, manual = TRUE, inherit.aes = FALSE,
              textsize = 3.4, tip_length = 0.01, vjust = -0.3,
              aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y_position)) +
  geom_text(data = ncount_b, inherit.aes = FALSE,
            aes(x = group, y = y, label = paste0("n=", n)),
            size = 3.2, vjust = 1, color = "grey25") +
  coord_cartesian(ylim = lay$ylim) +
  facet_wrap(~ stratum) +
  labs(x = NULL, y = "Mean expression (z-scored within dataset)",
       title = "Expression by module group (bulk), pooled over all cell types",
       subtitle = "Cell-type-agnostic groups (Shared = (mod, ct, dir) in both AD & SCZ in any cell type)",
       caption = "Whiskers capped at 97th pct for display. Brackets: Wilcoxon BH-adjusted p (*** <0.001, ** <0.01, * <0.05).") +
  theme_bw(base_size = 17) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16), strip.text = element_text(size = 16),
        plot.title = element_text(size = 19, face = "bold"), plot.subtitle = element_text(size = 15),
        legend.position = "none")

ggsave(file.path(fig_dir, "mod_group_expression_allct_bulk.pdf"), p_bulk, width = 10, height = 8)
ggsave(file.path(fig_dir, "mod_group_expression_allct_bulk.png"), p_bulk, width = 10, height = 8, dpi = 200)

# ---------------------------------------------------------------------------
# Single-nucleus figure: 4 dodged series per group (in-plot brackets ambiguous),
# with a companion Cliff's-delta heatmap (mirrors mod_group_expression_figure.R).
# ---------------------------------------------------------------------------
sn <- filter(pooled, dataset != "Bulk")
ylim_use <- quantile(sn$z, c(0.01, 0.98), na.rm = TRUE)
ncount_s <- sn %>% group_by(stratum, group) %>% summarise(n = n_distinct(gene), .groups = "drop")
rng <- ylim_use[2] - ylim_use[1]
ncount_s$y <- ylim_use[1] - 0.05 * rng
ylim_use[1] <- ylim_use[1] - 0.11 * rng

p_sn <- ggplot(sn, aes(x = group, y = z, fill = series)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.3) +
  geom_text(data = ncount_s, inherit.aes = FALSE,
            aes(x = group, y = y, label = paste0("n=", n)),
            size = 3.2, vjust = 1, color = "grey25") +
  coord_cartesian(ylim = ylim_use) +
  facet_wrap(~ stratum) +
  labs(x = NULL, y = "Mean expression (z-scored within dataset)",
       fill = "Dataset / condition",
       title = "Expression by module group (single-nucleus), pooled over all cell types",
       subtitle = "Cell-type-agnostic groups; every (gene x cell type) value is a point",
       caption = "See companion heatmap for Wilcoxon effect sizes (Cliff's delta + BH-adj p stars).") +
  theme_bw(base_size = 17) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16), strip.text = element_text(size = 16),
        plot.title = element_text(size = 19, face = "bold"), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "mod_group_expression_allct_sn.pdf"), p_sn, width = 11, height = 8)
ggsave(file.path(fig_dir, "mod_group_expression_allct_sn.png"), p_sn, width = 11, height = 8, dpi = 200)

p_sn_stats <- delta_heatmap(
  filter(wilcox_res, dataset != "Bulk"), "stratum",
  comp_levels = c("Shared vs Neither", "Unshared vs Neither", "Shared vs Unshared"),
  title = "Single-nucleus Wilcoxon effect sizes (companion to sn boxplot, all cell types pooled)",
  caption = "Tile = Cliff's delta; text = BH-adj p stars (*** <0.001, ** <0.01, * <0.05).")
ggsave(file.path(fig_dir, "mod_group_expression_allct_sn_stats.pdf"), p_sn_stats, width = 14, height = 4)
ggsave(file.path(fig_dir, "mod_group_expression_allct_sn_stats.png"), p_sn_stats, width = 14, height = 4, dpi = 200)
cat("done\n")

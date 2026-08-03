# Per-shared-MODULE figure for mean expression, COLLAPSED over cell type. Mirrors
# shared_module_featured_figure.R but with one panel per shared MODULE (12 unique
# modules) instead of one per (mod, cell type, direction) triplet. Each panel
# features that module's genes ("Module"), compared against ALL Unshared and ALL
# Neither genes (cell-type-agnostic groups; see group_overall in ../mod_group_defs.R).
# Single-nucleus expression is POOLED over all 23 cell types; bulk is one value per
# gene. Bulk and single-nucleus rendered as separate figures (own y-axis).
# Unshared/Neither are the same set in every panel (only the featured module changes).

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression/mean_expression_defs.R"))

fig_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression")

# Base: every gene's z-scored expression (sc pooled over cell types + bulk once)
# with its cell-type-agnostic group.
base <- bind_rows(
  sc_z %>% inner_join(gene_tbl, by = "gene") %>% mutate(series = paste(dataset, condition)),
  bulk_z %>% inner_join(gene_tbl, by = "gene") %>%
    mutate(dataset = "Bulk", condition = "Bulk", series = "Bulk")
) %>% mutate(group_ov = group_overall(mod))

shared_mods <- sort(unique(shared_keys$mod))   # 12 unique shared modules
panel_levels <- sprintf("mod %d", sort(match(shared_mods, these_mods)))

feat_long <- map_dfr(shared_mods, function(m) {
  m_disp <- match(m, these_mods)               # display id on the 1016-universe scale
  base %>%
    mutate(plot_group = case_when(
             mod == m            ~ "Module",
             group_ov == "Unshared" ~ "Unshared",
             group_ov == "Neither"  ~ "Neither",
             TRUE                ~ NA_character_)) %>%
    filter(!is.na(plot_group)) %>%
    mutate(panel = sprintf("mod %d", m_disp), feat_mod = m_disp)
})
feat_long$plot_group <- factor(feat_long$plot_group, levels = c("Module", "Unshared", "Neither"))
feat_long$panel <- factor(feat_long$panel, levels = panel_levels)

# Compact backing data: bulk long (small) + per-(panel, group, series) summary.
# The full single-nucleus long table (~19M rows) is not written -- it is the
# all-cell-type pooled expression (mod_group_expression_allct_long_sn.csv)
# re-grouped per featured module, so it is deterministic and large.
fwrite(filter(feat_long, dataset == "Bulk"),
       file.path(fig_dir, "shared_module_featured_allct_long_bulk.csv"))
feat_summary <- feat_long %>%
  group_by(panel, feat_mod, series, plot_group) %>%
  summarise(n_genes = n_distinct(gene), n_points = n(),
            median_z = median(z), mean_z = mean(z), .groups = "drop")
fwrite(feat_summary, file.path(fig_dir, "shared_module_featured_allct_summary.csv"))

# ---------------------------------------------------------------------------
# Wilcoxon per featured module x series, with Cliff's delta.
# ---------------------------------------------------------------------------
pairs <- list(c("Module", "Neither"), c("Module", "Unshared"), c("Unshared", "Neither"))

wilcox_one <- function(d, g1, g2) {
  x <- d$expr_scale[d$plot_group == g1]; y <- d$expr_scale[d$plot_group == g2]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  w <- wilcox.test(x, y); U <- unname(w$statistic)
  np <- as.numeric(length(x)) * as.numeric(length(y))   # numeric: avoid integer overflow at large n
  tibble(group1 = g1, group2 = g2, n1 = length(x), n2 = length(y),
         median1 = median(x), median2 = median(y), median_diff = median(x) - median(y),
         W = U, cliffs_delta = 2 * U / np - 1, p = w$p.value)
}

wilcox_res <- feat_long %>%
  group_by(panel, feat_mod, dataset, condition) %>%
  group_modify(~ map_dfr(pairs, function(pr) wilcox_one(.x, pr[1], pr[2]))) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH")) %>%
  arrange(feat_mod, dataset, condition, group1, group2)
fwrite(filter(wilcox_res, dataset != "Bulk"), file.path(fig_dir, "shared_module_featured_allct_wilcoxon_sn.csv"))
fwrite(filter(wilcox_res, dataset == "Bulk"), file.path(fig_dir, "shared_module_featured_allct_wilcoxon_bulk.csv"))

cat("panels:", nlevels(feat_long$panel), "| tests:", nrow(wilcox_res),
    "| BH<0.05:", sum(wilcox_res$p_adj_BH < 0.05), "\n")
md <- subset(wilcox_res, group1 == "Module" & group2 == "Neither")
cat("Module>Neither (positive Cliff's delta):", sum(md$cliffs_delta > 0), "/", nrow(md), "\n")

# ---------------------------------------------------------------------------
# Bulk figure: one panel per module, one box per group, Wilcoxon brackets.
# ---------------------------------------------------------------------------
bulk <- filter(feat_long, dataset == "Bulk")
cap  <- quantile(bulk$z, c(0.01, 0.97), na.rm = TRUE)
bulk <- mutate(bulk, z = pmin(pmax(z, cap[1]), cap[2]))
lay  <- bracket_layout(signif_brackets(filter(wilcox_res, dataset == "Bulk"),
                                       levels(feat_long$plot_group)),
                       bulk, "panel", "z", group_col = "plot_group")
ncount_b <- bulk %>% group_by(panel, plot_group) %>% summarise(n = n_distinct(gene), .groups = "drop")
rng <- lay$ylim[2] - lay$ylim[1]
ncount_b$y <- lay$ylim[1] - 0.05 * rng
lay$ylim[1] <- lay$ylim[1] - 0.11 * rng

p_bulk <- ggplot(bulk, aes(x = plot_group, y = z, fill = series)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) +
  geom_signif(data = lay$sig, manual = TRUE, inherit.aes = FALSE,
              textsize = 2.9, tip_length = 0.01, vjust = -0.3,
              aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y_position)) +
  geom_text(data = ncount_b, inherit.aes = FALSE,
            aes(x = plot_group, y = y, label = paste0("n=", n)),
            size = 2.7, vjust = 1, color = "grey25") +
  coord_cartesian(ylim = lay$ylim) +
  facet_wrap(~ panel, ncol = 3) +
  labs(x = NULL, y = "Mean expression (z-scored within dataset)",
       title = "Each shared module featured vs Unshared & Neither genes (bulk), all cell types pooled",
       subtitle = "'Module' = the featured shared module's genes; cell-type-agnostic groups",
       caption = "Whiskers capped at 97th pct for display. Brackets: Wilcoxon BH-adjusted p (*** <0.001, ** <0.01, * <0.05).") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15), strip.text = element_text(size = 13),
        plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 15),
        legend.position = "none")

ggsave(file.path(fig_dir, "shared_module_featured_allct_bulk.pdf"), p_bulk, width = 13, height = 16)
ggsave(file.path(fig_dir, "shared_module_featured_allct_bulk.png"), p_bulk, width = 13, height = 16, dpi = 200)

# ---------------------------------------------------------------------------
# Single-nucleus figure: 4 dodged series per group + companion Cliff's-delta heatmap.
# ---------------------------------------------------------------------------
sn <- filter(feat_long, dataset != "Bulk")
ylim_use <- quantile(sn$z, c(0.01, 0.98), na.rm = TRUE)
ncount_s <- sn %>% group_by(panel, plot_group) %>% summarise(n = n_distinct(gene), .groups = "drop")
rng <- ylim_use[2] - ylim_use[1]
ncount_s$y <- ylim_use[1] - 0.05 * rng
ylim_use[1] <- ylim_use[1] - 0.11 * rng

p_sn <- ggplot(sn, aes(x = plot_group, y = z, fill = series)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) +
  geom_text(data = ncount_s, inherit.aes = FALSE,
            aes(x = plot_group, y = y, label = paste0("n=", n)),
            size = 2.7, vjust = 1, color = "grey25") +
  coord_cartesian(ylim = ylim_use) +
  facet_wrap(~ panel, ncol = 3) +
  labs(x = NULL, y = "Mean expression (z-scored within dataset)",
       fill = "Dataset / condition",
       title = "Each shared module featured vs Unshared & Neither genes (single-nucleus), all cell types pooled",
       subtitle = "'Module' = the featured shared module's genes; every (gene x cell type) value is a point",
       caption = "See companion heatmap for Wilcoxon effect sizes (Cliff's delta + BH-adj p stars).") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15), strip.text = element_text(size = 13),
        plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 14), legend.text = element_text(size = 13),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "shared_module_featured_allct_sn.pdf"), p_sn, width = 13, height = 16)
ggsave(file.path(fig_dir, "shared_module_featured_allct_sn.png"), p_sn, width = 13, height = 16, dpi = 200)

sn_stats <- filter(wilcox_res, dataset != "Bulk")
p_sn_stats <- delta_heatmap(
  sn_stats, "panel", ylevels = rev(panel_levels),
  comp_levels = c("Module vs Neither", "Module vs Unshared", "Unshared vs Neither"),
  title = "Single-nucleus Wilcoxon effect sizes per featured module (companion to sn boxplots, all cell types pooled)",
  caption = "Tile = Cliff's delta; text = BH-adj p stars (*** <0.001, ** <0.01, * <0.05).")
ggsave(file.path(fig_dir, "shared_module_featured_allct_sn_stats.pdf"), p_sn_stats, width = 12, height = 7)
ggsave(file.path(fig_dir, "shared_module_featured_allct_sn_stats.png"), p_sn_stats, width = 12, height = 7, dpi = 200)
cat("done\n")

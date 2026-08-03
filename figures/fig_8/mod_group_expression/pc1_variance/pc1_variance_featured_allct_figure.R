# Per-shared-MODULE figure for PC1 variance, COLLAPSED over cell type. Mirrors
# pc1_variance_featured_figure.R but with one panel per shared MODULE (12 unique
# modules) instead of one per (mod, cell type, direction) triplet. Each panel
# features that module's genes ("Module"), with per-gene PC1 R² (bulk), compared
# against ALL Unshared and ALL Neither genes (cell-type-agnostic groups; see
# group_overall in ../mod_group_defs.R). Unshared/Neither are the same set in
# every panel (only the featured module changes).

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance/pc1_variance_defs.R"))

fig_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance")

# Per-gene R² (from pc1_variance_defs.R) + each gene's cell-type-agnostic group.
base <- r2_tbl %>% mutate(group_ov = group_overall(mod), series = "Bulk")

shared_mods <- sort(unique(shared_keys$mod))   # 12 unique shared modules
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
feat_long$panel <- factor(feat_long$panel,
                          levels = sprintf("mod %d", sort(match(shared_mods, these_mods))))
fwrite(feat_long, file.path(fig_dir, "pc1_variance_featured_allct_long.csv"))

# ---------------------------------------------------------------------------
# Wilcoxon per featured module, with Cliff's delta.
# ---------------------------------------------------------------------------
pairs <- list(c("Module", "Neither"), c("Module", "Unshared"), c("Unshared", "Neither"))

wilcox_one <- function(d, g1, g2) {
  x <- d$r2[d$plot_group == g1]; y <- d$r2[d$plot_group == g2]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  w <- wilcox.test(x, y); U <- unname(w$statistic)
  tibble(group1 = g1, group2 = g2, n1 = length(x), n2 = length(y),
         median1 = median(x), median2 = median(y), median_diff = median(x) - median(y),
         W = U, cliffs_delta = 2 * U / (length(x) * length(y)) - 1, p = w$p.value)
}

wilcox_res <- feat_long %>%
  mutate(dataset = "Bulk", condition = "Bulk") %>%
  group_by(panel, feat_mod, dataset, condition) %>%
  group_modify(~ map_dfr(pairs, function(pr) wilcox_one(.x, pr[1], pr[2]))) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH")) %>%
  arrange(feat_mod, group1, group2)
fwrite(wilcox_res, file.path(fig_dir, "pc1_variance_featured_allct_wilcoxon.csv"))

cat("panels:", nlevels(feat_long$panel), "| tests:", nrow(wilcox_res),
    "| BH<0.05:", sum(wilcox_res$p_adj_BH < 0.05), "\n")
md <- subset(wilcox_res, group1 == "Module" & group2 == "Neither")
cat("Module>Neither (positive Cliff's delta):", sum(md$cliffs_delta > 0), "/", nrow(md), "\n")

# ---------------------------------------------------------------------------
# Figure: one panel per shared module, bulk; Wilcoxon brackets above whiskers.
# ---------------------------------------------------------------------------
sig <- signif_brackets(wilcox_res, levels(feat_long$plot_group))
lay <- bracket_layout(sig, feat_long, "panel", "r2", group_col = "plot_group")

ncount <- feat_long %>% group_by(panel, plot_group) %>%
  summarise(n = n_distinct(gene), .groups = "drop")
rng <- lay$ylim[2] - lay$ylim[1]
ncount$y <- lay$ylim[1] - 0.05 * rng
lay$ylim[1] <- lay$ylim[1] - 0.11 * rng

p <- ggplot(feat_long, aes(x = plot_group, y = r2, fill = series)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) +
  geom_signif(data = lay$sig, manual = TRUE, inherit.aes = FALSE,
              textsize = 2.9, tip_length = 0.01, vjust = -0.3,
              aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y_position)) +
  geom_text(data = ncount, inherit.aes = FALSE,
            aes(x = plot_group, y = y, label = paste0("n=", n)),
            size = 2.7, vjust = 1, color = "grey25") +
  coord_cartesian(ylim = lay$ylim) +
  facet_wrap(~ panel, ncol = 3) +
  labs(x = NULL, y = "PC1 variance explained (%, per-gene R²)",
       fill = "Dataset",
       title = "Each shared module featured: per-gene PC1 variance explained vs Unshared & Neither genes",
       subtitle = "'Module' = the featured shared module's genes (bulk); cell-type-agnostic groups, all cell types pooled",
       caption = "Brackets: Wilcoxon BH-adjusted p (*** <0.001, ** <0.01, * <0.05).") +
  theme_bw(base_size = 15) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title   = element_text(size = 15),
        strip.text   = element_text(size = 13),
        plot.title   = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 15),
        legend.position = "none")

ggsave(file.path(fig_dir, "pc1_variance_featured_allct.pdf"), p, width = 13, height = 16)
ggsave(file.path(fig_dir, "pc1_variance_featured_allct.png"), p, width = 13, height = 16, dpi = 200)
cat("done\n")

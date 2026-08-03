# Per-shared-module figure for PC1 variance: feature each shared (module, cell
# type, direction) triplet, with that module's genes ("Module") compared against
# all Unshared and all Neither genes in that cell type. The value is each gene's
# % variance explained by its OWN module's PC1 (bulk; see pc1_variance_defs.R).
# Mirrors mean_expression/shared_module_featured_bulk but with per-gene PC1 R².
# One panel per shared triplet (18 total), bulk only.

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance/pc1_variance_defs.R"))

fig_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/pc1_variance")

dir_lab <- function(d) ifelse(d < 0, "down", "up")

# For each shared triplet, relabel that module's genes as "Module" and keep the
# Unshared / Neither genes of the same cell type.
feat_long <- pmap_dfr(shared_keys, function(mod, ct, Direction) {
  m <- mod; c <- ct; d <- Direction
  m_disp <- match(m, these_mods)            # display module id on the 1016-universe scale (these_mods sorted ascending)
  build_grouped_r2(c) %>%
    mutate(plot_group = case_when(
             mod == m            ~ "Module",
             group == "Unshared" ~ "Unshared",
             group == "Neither"  ~ "Neither",
             TRUE                ~ NA_character_)) %>%
    filter(!is.na(plot_group)) %>%
    mutate(panel = sprintf("mod %d · %s (%s)", m_disp, c, dir_lab(d)),
           feat_mod = m_disp, feat_ct = c)
})
feat_long$plot_group <- factor(feat_long$plot_group,
                               levels = c("Module", "Unshared", "Neither"))
fwrite(feat_long, file.path(fig_dir, "pc1_variance_featured_long.csv"))

# ---------------------------------------------------------------------------
# Wilcoxon tests per featured triplet, with Cliff's delta. Run before the figure
# so it can annotate the brackets.
# ---------------------------------------------------------------------------
pairs <- list(c("Module", "Neither"),
              c("Module", "Unshared"),
              c("Unshared", "Neither"))

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
  group_by(panel, feat_mod, feat_ct, dataset, condition) %>%
  group_modify(~ map_dfr(pairs, function(pr) wilcox_one(.x, pr[1], pr[2]))) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH")) %>%
  arrange(feat_mod, feat_ct, group1, group2)
fwrite(wilcox_res, file.path(fig_dir, "pc1_variance_featured_wilcoxon.csv"))

cat("panels:", length(unique(feat_long$panel)),
    "| tests:", nrow(wilcox_res),
    "| BH<0.05:", sum(wilcox_res$p_adj_BH < 0.05), "\n")
md <- subset(wilcox_res, group1 == "Module" & group2 == "Neither")
cat("Module>Neither (positive Cliff's delta):", sum(md$cliffs_delta > 0), "/", nrow(md), "\n")

# ---------------------------------------------------------------------------
# Figure: one panel per shared triplet, bulk; Wilcoxon brackets above the
# whiskers (single box per group; see bracket_layout in mod_group_defs.R).
# ---------------------------------------------------------------------------
sig <- signif_brackets(wilcox_res, levels(feat_long$plot_group))
lay <- bracket_layout(sig, feat_long, "panel", "r2", group_col = "plot_group")

# Genes analyzed per group, labelled in a reserved strip just below the boxes.
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
       title = "Each shared (module, cell type, direction) featured: per-gene PC1 variance explained vs Unshared & Neither genes",
       subtitle = "'Module' = the featured shared module's genes (bulk)",
       caption = "Brackets: Wilcoxon BH-adjusted p (*** <0.001, ** <0.01, * <0.05).") +
  theme_bw(base_size = 15) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title   = element_text(size = 15),
        strip.text   = element_text(size = 13),
        plot.title   = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 13),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "pc1_variance_featured.pdf"), p, width = 13, height = 21)
ggsave(file.path(fig_dir, "pc1_variance_featured.png"), p, width = 13, height = 21, dpi = 200)
cat("done\n")

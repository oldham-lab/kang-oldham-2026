# Per-shared-module figure: feature each shared (module, cell type, direction)
# triplet separately. "Shared" is replaced by the single featured module's
# genes ("Module"), compared against all Unshared and all Neither genes in that
# triplet's cell type (cell-type-resolved grouping; see mod_group_defs.R).
# One panel per shared triplet (18 total), each shown only in its own cell type.

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression/mean_expression_defs.R"))

fig_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression")

dir_lab <- function(d) ifelse(d < 0, "down", "up")

# For each shared triplet, relabel that module's genes as "Module" and keep the
# Unshared / Neither genes of the same cell type.
feat_long <- pmap_dfr(shared_keys, function(mod, ct, Direction) {
  m <- mod; c <- ct; d <- Direction
  m_disp <- match(m, these_mods)            # display module id on the 1016-universe scale (these_mods sorted ascending)
  build_grouped(c) %>%
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

fwrite(filter(feat_long, dataset != "Bulk"), file.path(fig_dir, "shared_module_featured_long_sn.csv"))
fwrite(filter(feat_long, dataset == "Bulk"), file.path(fig_dir, "shared_module_featured_long_bulk.csv"))

# ---------------------------------------------------------------------------
# Wilcoxon tests per featured triplet x dataset/condition, with Cliff's delta.
# Run before the figures so the bulk figure can annotate these results.
# ---------------------------------------------------------------------------
pairs <- list(c("Module", "Neither"),
              c("Module", "Unshared"),
              c("Unshared", "Neither"))

wilcox_one <- function(d, g1, g2) {
  x <- d$expr_scale[d$plot_group == g1]
  y <- d$expr_scale[d$plot_group == g2]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  w <- wilcox.test(x, y)
  U <- unname(w$statistic)
  tibble(group1 = g1, group2 = g2,
         n1 = length(x), n2 = length(y),
         median1 = median(x), median2 = median(y),
         median_diff = median(x) - median(y),
         W = U, cliffs_delta = 2 * U / (length(x) * length(y)) - 1,
         p = w$p.value)
}

wilcox_res <- feat_long %>%
  group_by(panel, feat_mod, feat_ct, dataset, condition) %>%
  group_modify(~ map_dfr(pairs, function(pr) wilcox_one(.x, pr[1], pr[2]))) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH")) %>%
  arrange(feat_mod, feat_ct, dataset, condition, group1, group2)

fwrite(filter(wilcox_res, dataset != "Bulk"), file.path(fig_dir, "shared_module_featured_wilcoxon_sn.csv"))
fwrite(filter(wilcox_res, dataset == "Bulk"), file.path(fig_dir, "shared_module_featured_wilcoxon_bulk.csv"))

cat("panels:", length(unique(feat_long$panel)),
    "| tests:", nrow(wilcox_res),
    "| BH<0.05:", sum(wilcox_res$p_adj_BH < 0.05), "\n")
md <- subset(wilcox_res, group1 == "Module" & group2 == "Neither")
cat("Module>Neither (positive Cliff's delta):", sum(md$cliffs_delta > 0), "/", nrow(md), "\n")

# ---------------------------------------------------------------------------
# Figures: one panel per shared triplet, bulk and single-nucleus rendered as
# separate figures (each with its own y-axis zoom). `sig` (optional) draws
# pairwise Wilcoxon brackets (BH-adj p + stars); used for bulk (single box/group).
# ---------------------------------------------------------------------------
make_fig <- function(d, subtitle, sig = NULL) {
  sig_layer <- NULL
  if (!is.null(sig) && nrow(sig) > 0) {
    cap <- quantile(d$z, c(0.01, 0.97), na.rm = TRUE)        # winsorize whiskers for
    d   <- mutate(d, z = pmin(pmax(z, cap[1]), cap[2]))      # display so boxes stay legible
    lay <- bracket_layout(sig, d, "panel", "z", group_col = "plot_group")
    sig_layer <- geom_signif(data = lay$sig, manual = TRUE, inherit.aes = FALSE,
                             textsize = 2.9, tip_length = 0.01, vjust = -0.3,
                             aes(xmin = xmin, xmax = xmax,
                                 annotations = label, y_position = y_position))
    ylim_use <- lay$ylim
  } else {
    ylim_use <- quantile(d$z, c(0.01, 0.98), na.rm = TRUE)   # sn: zoom unchanged
  }
  # Genes analyzed per group (identical across series within a figure), labelled
  # in a reserved strip just below the boxes in each panel.
  ncount <- d %>% group_by(panel, plot_group) %>%
    summarise(n = n_distinct(gene), .groups = "drop")
  rng <- ylim_use[2] - ylim_use[1]
  ncount$y <- ylim_use[1] - 0.05 * rng
  ylim_use[1] <- ylim_use[1] - 0.11 * rng
  count_layer <- geom_text(data = ncount, inherit.aes = FALSE,
                           aes(x = plot_group, y = y, label = paste0("n=", n)),
                           size = 2.7, vjust = 1, color = "grey25")
  ggplot(d, aes(x = plot_group, y = z, fill = series)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.25) +
    sig_layer +
    count_layer +
    coord_cartesian(ylim = ylim_use) +
    facet_wrap(~ panel, ncol = 3) +
    labs(x = NULL, y = "Mean expression (z-scored within dataset)",
         fill = "Dataset / condition",
         title = "Each shared (module, cell type, direction) featured vs Unshared & Neither genes in that cell type",
         subtitle = subtitle,
         caption = "'Module' = the featured shared module's genes. Outliers hidden; whiskers capped at 97th pct for display. Brackets: Wilcoxon BH-adjusted p (*** <0.001, ** <0.01, * <0.05).") +
    theme_bw(base_size = 15) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 13),
          axis.text.y  = element_text(size = 13),
          axis.title   = element_text(size = 15),
          strip.text   = element_text(size = 13),
          plot.title   = element_text(size = 17, face = "bold"),
          plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 13),
          legend.position = "bottom")
}

sn_feat   <- filter(feat_long, dataset != "Bulk")
bulk_feat <- filter(feat_long, dataset == "Bulk")

bulk_sig <- signif_brackets(filter(wilcox_res, dataset == "Bulk"),
                            levels(feat_long$plot_group))

p_sn   <- make_fig(sn_feat,   "Single-nucleus (Liu, Gabitto; Con & allAD)")
p_bulk <- make_fig(bulk_feat, "Bulk (all samples)", sig = bulk_sig)

ggsave(file.path(fig_dir, "shared_module_featured_sn.pdf"),   p_sn,   width = 13, height = 21)
ggsave(file.path(fig_dir, "shared_module_featured_sn.png"),   p_sn,   width = 13, height = 21, dpi = 200)
ggsave(file.path(fig_dir, "shared_module_featured_bulk.pdf"), p_bulk, width = 13, height = 21)
ggsave(file.path(fig_dir, "shared_module_featured_bulk.png"), p_bulk, width = 13, height = 21, dpi = 200)

# Companion effect-size heatmap for the single-nucleus per-triplet figure
# (rows = featured triplet; see delta_heatmap in mod_group_defs.R).
sn_stats <- filter(wilcox_res, dataset != "Bulk")
p_sn_stats <- delta_heatmap(
  sn_stats, "panel",
  ylevels = rev(unique(sn_stats$panel[order(sn_stats$feat_mod, sn_stats$feat_ct)])),
  comp_levels = c("Module vs Neither", "Module vs Unshared", "Unshared vs Neither"),
  title = "Single-nucleus Wilcoxon effect sizes per featured triplet (companion to sn boxplots)",
  caption = "Tile = Cliff's delta; text = BH-adj p stars (*** <0.001, ** <0.01, * <0.05).")
ggsave(file.path(fig_dir, "shared_module_featured_sn_stats.pdf"), p_sn_stats, width = 12, height = 9)
ggsave(file.path(fig_dir, "shared_module_featured_sn_stats.png"), p_sn_stats, width = 12, height = 9, dpi = 200)
cat("done\n")

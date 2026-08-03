# Figure + Wilcoxon: gene expression across three module groups, where
# "Shared" is defined per (module, cell type, direction) -- see mod_group_defs.R.
# Grouping is cell-type-resolved and shown for ALL cell types; cell types with
# no shared triplet show only Unshared / Neither (no Shared box).
#
# Bulk and single-nucleus results are rendered as SEPARATE figures (each with its
# own y-axis zoom): single-nucleus = Liu (MIT) + Gabitto (SeAAD2024), Con & allAD;
# bulk = all samples (celltype-agnostic, re-grouped by each cell type's membership).

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression/mean_expression_defs.R"))

fig_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mean_expression")

# build_grouped() and the cell-type-resolved grouping live in mod_group_defs.R
all_long <- build_grouped(expr_celltypes)
fwrite(filter(all_long, dataset != "Bulk"), file.path(fig_dir, "mod_group_expression_long_sn.csv"))
fwrite(filter(all_long, dataset == "Bulk"), file.path(fig_dir, "mod_group_expression_long_bulk.csv"))
cat("rows per group:\n"); print(table(all_long$group))

# ---------------------------------------------------------------------------
# Wilcoxon rank-sum tests between groups (per dataset x condition x cell type)
# with Cliff's delta effect size (= rank-biserial r, range [-1, 1]).
# Run before the figures so the bulk figure can annotate these results.
# ---------------------------------------------------------------------------
pairs <- list(c("Shared", "Neither"),
              c("Unshared", "Neither"),
              c("Shared", "Unshared"))

wilcox_one <- function(d, g1, g2) {
  x <- d$expr_scale[d$group == g1]
  y <- d$expr_scale[d$group == g2]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  w <- wilcox.test(x, y)
  U <- unname(w$statistic)                              # Mann-Whitney U for x (group1)
  cliffs_delta <- 2 * U / (length(x) * length(y)) - 1   # = rank-biserial r, [-1, 1]
  tibble(group1 = g1, group2 = g2,
         n1 = length(x), n2 = length(y),
         median1 = median(x), median2 = median(y),
         median_diff = median(x) - median(y),
         W = U, cliffs_delta = cliffs_delta, p = w$p.value)
}

wilcox_res <- all_long %>%       # build_grouped keeps expr_scale from sc_z / bulk_z
  group_by(dataset, condition, celltype) %>%
  group_modify(~ map_dfr(pairs, function(pr) wilcox_one(.x, pr[1], pr[2]))) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH")) %>%
  arrange(dataset, condition, celltype, group1, group2)

fwrite(filter(wilcox_res, dataset != "Bulk"), file.path(fig_dir, "mod_group_expression_wilcoxon_sn.csv"))
fwrite(filter(wilcox_res, dataset == "Bulk"), file.path(fig_dir, "mod_group_expression_wilcoxon_bulk.csv"))

cat("\nWilcoxon tests:", nrow(wilcox_res),
    "| BH<0.05:", sum(wilcox_res$p_adj_BH < 0.05), "\n")
cat("Shared>Neither (positive Cliff's delta):",
    with(subset(wilcox_res, group1=="Shared" & group2=="Neither"), sum(cliffs_delta > 0)),
    "/", nrow(subset(wilcox_res, group1=="Shared" & group2=="Neither")), "\n")

# ---------------------------------------------------------------------------
# Figures: bulk and single-nucleus rendered separately, each with its own
# y-axis zoom (all cell types faceted; non-shared cell types lack a Shared box).
# `sig` (optional) draws pairwise Wilcoxon brackets (BH-adj p + stars); used for
# bulk, where each group has a single box so brackets are unambiguous.
# ---------------------------------------------------------------------------
make_fig <- function(d, subtitle, sig = NULL) {
  sig_layer <- NULL
  if (!is.null(sig) && nrow(sig) > 0) {
    cap <- quantile(d$z, c(0.01, 0.97), na.rm = TRUE)        # winsorize whiskers for
    d   <- mutate(d, z = pmin(pmax(z, cap[1]), cap[2]))      # display so boxes stay legible
    lay <- bracket_layout(sig, d, "celltype", "z")           # brackets above capped whiskers
    sig_layer <- geom_signif(data = lay$sig, manual = TRUE, inherit.aes = FALSE,
                             textsize = 3.0, tip_length = 0.01, vjust = -0.3,
                             aes(xmin = xmin, xmax = xmax,
                                 annotations = label, y_position = y_position))
    ylim_use <- lay$ylim
  } else {
    ylim_use <- quantile(d$z, c(0.01, 0.98), na.rm = TRUE)   # sn: zoom unchanged
  }
  # Genes analyzed per group (identical across series within a figure), labelled
  # in a reserved strip just below the boxes in each facet.
  ncount <- d %>% group_by(celltype, group) %>%
    summarise(n = n_distinct(gene), .groups = "drop")
  rng <- ylim_use[2] - ylim_use[1]
  ncount$y <- ylim_use[1] - 0.05 * rng
  ylim_use[1] <- ylim_use[1] - 0.11 * rng
  count_layer <- geom_text(data = ncount, inherit.aes = FALSE,
                           aes(x = group, y = y, label = paste0("n=", n)),
                           size = 2.7, vjust = 1, color = "grey25")
  ggplot(d, aes(x = group, y = z, fill = series)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.3) +
    sig_layer +
    count_layer +
    coord_cartesian(ylim = ylim_use) +
    facet_wrap(~ celltype) +
    labs(x = NULL, y = "Mean expression (z-scored within dataset)",
         fill = "Dataset / condition",
         title = "Expression by module group, Shared = (module, cell type, direction) in both AD & SCZ",
         subtitle = subtitle,
         caption = "All cell types. Outliers hidden; whiskers capped at 97th pct for display. Brackets: Wilcoxon BH-adjusted p (*** <0.001, ** <0.01, * <0.05).") +
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
}

sn_long   <- filter(all_long, dataset != "Bulk")
bulk_long <- filter(all_long, dataset == "Bulk")

bulk_sig <- signif_brackets(filter(wilcox_res, dataset == "Bulk"), levels(all_long$group))

p_sn   <- make_fig(sn_long,   "Single-nucleus (Liu, Gabitto; Con & allAD)")
p_bulk <- make_fig(bulk_long, "Bulk (all samples; celltype-agnostic, re-grouped per cell type)", sig = bulk_sig)

ggsave(file.path(fig_dir, "mod_group_expression_by_celltype_sn.pdf"),   p_sn,   width = 16, height = 16)
ggsave(file.path(fig_dir, "mod_group_expression_by_celltype_sn.png"),   p_sn,   width = 16, height = 16, dpi = 200)
ggsave(file.path(fig_dir, "mod_group_expression_by_celltype_bulk.pdf"), p_bulk, width = 16, height = 16)
ggsave(file.path(fig_dir, "mod_group_expression_by_celltype_bulk.png"), p_bulk, width = 16, height = 16, dpi = 200)

# Companion effect-size heatmap for the single-nucleus figure (4 dodged series
# per group make in-plot brackets ambiguous; see delta_heatmap in mod_group_defs.R).
p_sn_stats <- delta_heatmap(
  filter(wilcox_res, dataset != "Bulk"), "celltype", ylevels = rev(expr_celltypes),
  comp_levels = c("Shared vs Neither", "Unshared vs Neither", "Shared vs Unshared"),
  title = "Single-nucleus Wilcoxon effect sizes (companion to sn boxplots)",
  caption = "Tile = Cliff's delta; text = BH-adj p stars (*** <0.001, ** <0.01, * <0.05). Empty cell = group absent in that cell type.")
ggsave(file.path(fig_dir, "mod_group_expression_by_celltype_sn_stats.pdf"), p_sn_stats, width = 13, height = 9)
ggsave(file.path(fig_dir, "mod_group_expression_by_celltype_sn_stats.png"), p_sn_stats, width = 13, height = 9, dpi = 200)
cat("done\n")

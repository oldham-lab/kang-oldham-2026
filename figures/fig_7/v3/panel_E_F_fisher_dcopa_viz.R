# Per-celltype Fisher overlap p-value: small-multiple 8×8 heatmaps

library(data.table)
library(tidyverse)
library(showtext)
library(cowplot)
library(svglite)
showtext_auto()

save_plot <- function(path_no_ext, plot_expr, width, height, svg_bg = "white") {
  cairo_pdf(paste0(path_no_ext, ".pdf"), width = width, height = height)
  print(plot_expr)
  dev.off()
  svglite(paste0(path_no_ext, ".svg"), width = width, height = height, bg = svg_bg)
  print(plot_expr)
  dev.off()
  message("Saved: ", basename(path_no_ext), ".pdf / .svg")
}
showtext_auto()

# ============================================================
# PLACEHOLDERS
# ============================================================
FULL_TABLE_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3/panel_B_fisher_dcopa_full.csv")
DCOPA_PATH      <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3/panel_B_dcopa_genelist.csv")
SAVE_DIR        <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3")

FDR_THRESH  <- 0.05
NEGLOG_CAP  <- 10

# ============================================================
# SET METADATA
# ============================================================
SET_KEYS <- c(
  "Gabitto_AllADVsCon_DFC_ROSMAP",
  "Liu_AllADVsCon_DFC_ROSMAP",
  "Gabitto_AllADVsCon_MTG_ROSMAP",
  "Liu_AllADVsCon_MTG_ROSMAP",
  "Gabitto_AllADVsCon_DFC",
  "Liu_AllADVsCon_DFC",
  "Gabitto_AllADVsCon_MTG",
  "Liu_AllADVsCon_MTG"
)

abbrev <- c(
  "Gabitto_AllADVsCon_DFC_ROSMAP" = "Gab_DFC_ADmods",
  "Liu_AllADVsCon_DFC_ROSMAP"     = "Liu_DFC_ADmods",
  "Gabitto_AllADVsCon_MTG_ROSMAP" = "Gab_MTG_ADmods",
  "Liu_AllADVsCon_MTG_ROSMAP"     = "Liu_MTG_ADmods",
  "Gabitto_AllADVsCon_DFC"        = "Gab_DFC",
  "Liu_AllADVsCon_DFC"            = "Liu_DFC",
  "Gabitto_AllADVsCon_MTG"        = "Gab_MTG",
  "Liu_AllADVsCon_MTG"            = "Liu_MTG"
)

# ============================================================
# CELLTYPE MAPPING  (Liu MIT names → SEA canonical names)
# ============================================================
map_list <- list(
    c("Endothelial", "SMC", "VLMC", "End", "Per"),
    c("L4 IT", "Exc L4-5 IT-2", "Exc L3-4 IT","Exc L4-5 IT-1"),
    c("L5 ET",  "Exc L5 ET"),
    c("L5 IT", "Exc L4-5 IT-2", "Exc L4-5 IT-1","Exc L3-5 IT", "Exc L5-6 IT"),
    c("L5/6 NP", "Exc L5/6 NP"),
    c("Lamp5", "Inh LAMP5"),
    c("Pvalb", "Inh PVALB"),
    c("Sst", "Inh SST"),
    c("L6 IT", "Exc L5-6 IT", "Exc L6 IT"),
    c("L6 IT Car3", "Exc L5/6 IT Car3"),
    c("L6 CT", "Exc L6 CT"),
    c("Pax6","Inh PAX6"),
    c("Astrocyte", "Ast"),
    c("OPC", "OPC"),
    c("Vip",  "Inh VIP"),
    c("L6b", "Exc L6b"),
    c("L2/3 IT", "Exc L2-3 IT"),
    c("Microglia-PVM", "Mic"),
    c("Oligodendrocyte", "Oli")
    )

# Build order-independent canonical lookup after data is loaded (see below).
# For each synonym group, whichever member appears in the Gabitto data becomes
# the canonical name; all other members map to it.
build_canonical_map <- function(map_list, gabitto_celltypes) {
  lookup <- c()
  for (grp in map_list) {
    in_gab    <- grp[grp %in% gabitto_celltypes]
    canonical <- if (length(in_gab) >= 1) in_gab[1] else grp[1]
    lookup    <- c(lookup, setNames(rep(canonical, length(grp)), grp))
  }
  lookup
}

# ============================================================
# LOAD DATA
# ============================================================
message("Loading data ...")

ft <- fread(FULL_TABLE_PATH, data.table = FALSE) |>
  dplyr::mutate(
    neglog10_padj = pmin(-log10(p_adj + 1e-300), NEGLOG_CAP),
    sig           = p_adj < FDR_THRESH
  )

# ============================================================
# PREPARE PLOT DATA
# ============================================================

# Aggregate across directions: mean neglog10, any-sig flag
ft_ct <- ft |>
  dplyr::group_by(Group1, Group2, Celltype) |>
  dplyr::summarise(
    mean_neglog10 = mean(neglog10_padj, na.rm = TRUE),
    min_padj      = min(p_adj,          na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    sig   = min_padj < FDR_THRESH,
    ab_G1 = abbrev[Group1],
    ab_G2 = abbrev[Group2]
  )

# Drop celltypes with no tests involving a Liu group (Gabitto-only subtypes
# with no Liu equivalent)
liu_keys   <- SET_KEYS[grepl("Liu", SET_KEYS)]
ct_has_liu <- ft_ct |>
  dplyr::filter(Group1 %in% liu_keys | Group2 %in% liu_keys) |>
  dplyr::pull(Celltype) |>
  unique()
ft_ct <- dplyr::filter(ft_ct, Celltype %in% ct_has_liu)

# Celltype order: fixed display order (celltypes not in data are dropped)
ct_order_requested <- c(
  "L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT",
  "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6",
  "Pvalb", "Sncg", "Sst", "Vip", "Sst Chodl", "Astrocyte",
  "Oligodendrocyte", "OPC", "Microglia-PVM", "Endothelial", "VLMC"
)
ct_in_data <- unique(ft_ct$Celltype)
ct_order   <- c(
  intersect(ct_order_requested, ct_in_data),
  setdiff(ct_in_data, ct_order_requested)
)

# ── Direction-split data for Plot 1 ──────────────────────────────────────────
make_dir_ft <- function(ft_raw, direction_pattern) {
  ft_raw |>
    dplyr::filter(grepl(direction_pattern, Direction, ignore.case = TRUE)) |>
    dplyr::group_by(Group1, Group2, Celltype) |>
    dplyr::summarise(
      neglog10 = mean(neglog10_padj, na.rm = TRUE),
      min_padj = min(p_adj,          na.rm = TRUE),
      .groups  = "drop"
    ) |>
    dplyr::mutate(
      sig   = min_padj < FDR_THRESH,
      ab_G1 = abbrev[Group1],
      ab_G2 = abbrev[Group2]
    ) |>
    dplyr::filter(Celltype %in% ct_has_liu)
}

ft_higher <- make_dir_ft(ft, "higher")
ft_lower  <- make_dir_ft(ft, "lower")

make_hm_df <- function(ft_dir) {
  dplyr::bind_rows(
    ft_dir |> dplyr::rename(row = ab_G1, col = ab_G2),
    ft_dir |> dplyr::rename(row = ab_G2, col = ab_G1)
  ) |>
    dplyr::mutate(
      row      = factor(row,      levels = abbrev[SET_KEYS]),
      col      = factor(col,      levels = abbrev[SET_KEYS]),
      Celltype = factor(Celltype, levels = ct_order)
    )
}

hm_higher <- make_hm_df(ft_higher)
hm_lower  <- make_hm_df(ft_lower)

# Keep mean-aggregated hm_df for Plot 2
hm_df <- dplyr::bind_rows(
  ft_ct |> dplyr::rename(row = ab_G1, col = ab_G2),
  ft_ct |> dplyr::rename(row = ab_G2, col = ab_G1)
) |>
  dplyr::mutate(
    row      = factor(row,      levels = abbrev[SET_KEYS]),
    col      = factor(col,      levels = abbrev[SET_KEYS]),
    Celltype = factor(Celltype, levels = ct_order)
  )

# ============================================================
# PLOT 1: SMALL-MULTIPLE 8×8 HEATMAPS (higher & lower separate)
# ============================================================
message("Plotting small-multiple heatmaps ...")

n_cols <- 5
n_rows <- ceiling(length(ct_order) / n_cols)

base_theme <- theme_bw(base_size = 28) +
  theme(
    axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 22),
    axis.text.y   = element_text(size = 22),
    strip.text    = element_text(size = 24, face = "bold"),
    legend.title  = element_text(size = 22),
    legend.text   = element_text(size = 20),
    panel.spacing = unit(10, "pt")
  )

build_sm_plot <- function(hm, show_legend = FALSE) {
  sig_rows <- dplyr::filter(hm, sig)
  ggplot(hm, aes(x = col, y = row, fill = neglog10)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(
      data = sig_rows, aes(label = "*"),
      colour = "white", size = 12, vjust = 0.75, inherit.aes = TRUE
    ) +
    scale_fill_gradient(
      low = "#F7FBFF", high = "#08306B",
      limits   = c(0, NEGLOG_CAP),
      name     = "-log10(FDR)",
      na.value = "grey85"
    ) +
    facet_wrap(~ Celltype, ncol = n_cols) +
    coord_fixed() +
    base_theme +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title      = element_blank(),
      plot.subtitle   = element_blank()
    ) +
    labs(x = NULL, y = NULL)
}

if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)

pdf_w <- 22
pdf_h <- n_rows * 5.5 + 2.5

save_plot(file.path(SAVE_DIR, "panel_F"), build_sm_plot(hm_higher), pdf_w, pdf_h)
save_plot(file.path(SAVE_DIR, "panel_E"), build_sm_plot(hm_lower),  pdf_w, pdf_h)

# Legend as standalone PDF + SVG
p_for_legend <- build_sm_plot(hm_higher, show_legend = TRUE)
legend_grob  <- cowplot::get_legend(p_for_legend)
p_legend_out <- cowplot::ggdraw() + cowplot::draw_grob(legend_grob)
save_plot(file.path(SAVE_DIR, "panel_E_F_legend"), p_legend_out, 4, 3, svg_bg = "transparent")

# ============================================================
# PLOT 2: CELLTYPE × CELLTYPE HEATMAP PER DATASET PAIR
# ============================================================
message("Plotting celltype-by-celltype heatmaps per dataset pair ...")

# Deduplicate dataset pairs: keep canonical direction (Group1 <= Group2)
ct_per_pair <- ft_ct |>
  dplyr::filter(Group1 <= Group2) |>
  dplyr::mutate(pair_label = paste0(ab_G1, " vs ", ab_G2))

# Self-join to get all celltype × celltype combinations per dataset pair
ct_by_ct <- ct_per_pair |>
  dplyr::select(pair_label, Celltype, mean_neglog10, sig) |>
  dplyr::inner_join(
    ct_per_pair |>
      dplyr::select(pair_label, Celltype, mean_neglog10, sig) |>
      dplyr::rename(Celltype2 = Celltype, mean_neglog10_2 = mean_neglog10, sig2 = sig),
    by = "pair_label"
  ) |>
  dplyr::mutate(
    joint_neglog10 = (mean_neglog10 + mean_neglog10_2) / 2,
    joint_sig      = sig & sig2,
    ct1 = factor(Celltype,  levels = ct_order),
    ct2 = factor(Celltype2, levels = ct_order)
  )

sig_ct_df <- ct_by_ct |> dplyr::filter(joint_sig)

n_pairs     <- length(unique(ct_by_ct$pair_label))
n_pair_cols <- min(4, n_pairs)
n_pair_rows <- ceiling(n_pairs / n_pair_cols)
n_cts       <- length(ct_order)

p_ct <- ggplot(ct_by_ct, aes(x = ct2, y = ct1, fill = joint_neglog10)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(
    data       = sig_ct_df,
    aes(label  = "*"),
    colour     = "white", size = 6, vjust = 0.75, inherit.aes = TRUE
  ) +
  scale_fill_gradient(
    low      = "#F7FBFF",
    high     = "#08306B",
    limits   = c(0, NEGLOG_CAP),
    name     = "Mean\n-log10(FDR)",
    na.value = "grey85"
  ) +
  facet_wrap(~ pair_label, ncol = n_pair_cols) +
  coord_fixed() +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y   = element_text(size = 14),
    strip.text    = element_text(size = 16, face = "bold"),
    legend.title  = element_text(size = 16),
    legend.text   = element_text(size = 14),
    plot.title    = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 16),
    legend.position = "right",
    panel.spacing = unit(8, "pt")
  ) +
  labs(
    title    = "Pairwise Fisher overlap: celltype \u00d7 celltype per dataset pair",
    subtitle = paste0(
      "Asterisk = both celltypes FDR < ", FDR_THRESH, "; ",
      "fill = mean -log10(FDR) of the two celltypes; grey = not tested"
    ),
    x = NULL, y = NULL
  )

save_plot(file.path(SAVE_DIR, "panel_E_F_ct_by_ct"), p_ct,
          n_cts * 0.55 * n_pair_cols + 4,
          n_cts * 0.55 * n_pair_rows + 3)

library(data.table)
library(tidyverse)
library(patchwork)
library(cowplot)

# ============================================================
# PLACEHOLDERS
# ============================================================
SUMMARY_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "dcopa_overlap/dcopa_overlap_summary.csv")
GENES_PATH   <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "dcopa_overlap/dcopa_overlap_genes.csv")
SAVE_DIR     <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "dcopa_overlap/plots/")

FIG1_W <- 28
FIG1_H <- 12

FIG2_W <- 7
FIG2_H <- 9

# ============================================================
# FACTORS & COLOURS
# ============================================================
FACTOR_LEVELS <- c("MIT edgeR", "MIT DESeq2", "SEA edgeR", "SEA DESeq2", "dcopa")

FACTOR_COLOURS <- c(
  "MIT edgeR"  = "#4E79A7",
  "MIT DESeq2" = "#A0CBE8",
  "SEA edgeR"  = "#F28E2B",
  "SEA DESeq2" = "#FFBE7D",
  "dcopa"      = "#27AE60"
)

DATASETS <- c("MIT_edgeR", "MIT_DESeq2", "SEA_edgeR", "SEA_DESeq2")

# ============================================================
# LOAD DATA
# ============================================================
message("Loading data ...")
summary_raw <- fread(SUMMARY_PATH, data.table = FALSE)
genes_raw   <- fread(GENES_PATH,   data.table = FALSE)

if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)

all_cts <- sort(unique(c(summary_raw$celltype, genes_raw$celltype)))

# ============================================================
# BUILD LONG TABLES
# ============================================================
de_long <- summary_raw |>
  tidyr::pivot_longer(
    cols      = matches("^n_DE_(higher|lower)_"),
    names_to  = "key",
    values_to = "n_genes"
  ) |>
  dplyr::mutate(
    direction = ifelse(grepl("higher", key), "Higher in AD", "Lower in AD"),
    factor    = stringr::str_extract(key, paste(DATASETS, collapse = "|")),
    factor    = stringr::str_replace_all(factor, "_", " ")
  ) |>
  dplyr::select(celltype, factor, direction, n_genes)

dcopa_long <- summary_raw |>
  dplyr::select(celltype, n_dcopa_higher, n_dcopa_lower) |>
  tidyr::pivot_longer(
    cols      = c(n_dcopa_higher, n_dcopa_lower),
    names_to  = "direction",
    values_to = "n_genes"
  ) |>
  dplyr::mutate(
    direction = ifelse(direction == "n_dcopa_higher", "Higher in AD", "Lower in AD"),
    factor    = "dcopa"
  ) |>
  dplyr::select(celltype, factor, direction, n_genes)

plot_data <- dplyr::bind_rows(de_long, dcopa_long) |>
  dplyr::mutate(
    factor    = factor(factor,    levels = FACTOR_LEVELS),
    direction = factor(direction, levels = c("Higher in AD", "Lower in AD")),
    celltype  = factor(celltype,  levels = rev(all_cts))
  ) |>
  dplyr::filter(!is.na(n_genes))

pct_long <- summary_raw |>
  tidyr::pivot_longer(
    cols      = matches("^pct_overlap_(higher|lower)_"),
    names_to  = "key",
    values_to = "pct_overlap"
  ) |>
  dplyr::mutate(
    direction = ifelse(grepl("higher", key), "Higher in AD", "Lower in AD"),
    factor    = stringr::str_extract(key, paste(DATASETS, collapse = "|")),
    factor    = stringr::str_replace_all(factor, "_", " "),
    factor    = factor(factor,    levels = FACTOR_LEVELS),
    direction = factor(direction, levels = c("Higher in AD", "Lower in AD")),
    celltype  = factor(celltype,  levels = rev(all_cts))
  ) |>
  dplyr::select(celltype, factor, direction, pct_overlap) |>
  dplyr::filter(!is.na(pct_overlap))

# ============================================================
# SHARED THEME & SCALES
# ============================================================
base_theme <- theme_bw(base_size = 13) +
  theme(
    strip.background   = element_rect(fill = "grey94", colour = "grey70"),
    strip.text         = element_text(face = "bold", size = 13),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 12),
    axis.text.x        = element_text(size = 11),
    axis.title.x       = element_text(size = 12),
    legend.position    = "bottom",
    legend.key.size    = unit(0.5, "cm"),
    legend.text        = element_text(size = 12),
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 11, colour = "grey40")
  )

fill_scale <- scale_fill_manual(values = FACTOR_COLOURS, name = NULL, drop = FALSE)


# ============================================================
# FIGURE 1
# 8 columns: for each of the 4 DE datasets, a counts col + a % overlap col
# 2 rows: Higher in AD (top), Lower in AD (bottom)
# y-axis labels only on the leftmost column
# Column headers (dataset name) shown as plot titles on the top row only
# ============================================================
message("Building Figure 1: 8-column layout ...")

# Panel builders ---------------------------------------------------------------

# Counts panel: DE dataset bars + dcopa bars, 2 bars per celltype
make_counts_col <- function(ds, dir_label, show_y = FALSE, show_title = FALSE) {
  ds_label <- stringr::str_replace_all(ds, "_", " ")
  sub <- plot_data |>
    dplyr::filter(
      direction == dir_label,
      factor %in% c(ds_label, "dcopa")
    ) |>
    dplyr::mutate(factor = droplevels(factor))

  ggplot(sub, aes(x = n_genes, y = celltype, fill = factor)) +
    geom_col(position = position_dodge(width = 0.75),
             width = 0.7, colour = NA) +
    scale_fill_manual(
      values = FACTOR_COLOURS[c(ds_label, "dcopa")],
      name   = NULL, drop = FALSE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = if (show_title) ds_label else NULL,
      x     = "# genes",
      y     = NULL
    ) +
    base_theme +
    theme(
      plot.title      = element_text(face = "bold", size = 12,
                                     colour = FACTOR_COLOURS[[ds_label]]),
      axis.text.y     = if (show_y) element_text(size = 12) else element_blank(),
      axis.ticks.y    = if (show_y) element_line() else element_blank(),
      legend.position = "none"
    )
}

# % overlap panel: single bar per celltype for this dataset
make_pct_col <- function(ds, dir_label, show_y = FALSE) {
  ds_label <- stringr::str_replace_all(ds, "_", " ")
  sub <- pct_long |>
    dplyr::filter(direction == dir_label, factor == ds_label)

  ggplot(sub, aes(x = pct_overlap, y = celltype, fill = factor)) +
    geom_col(width = 0.5, colour = NA) +
    geom_text(
      aes(label = ifelse(pct_overlap > 0, sprintf("%.0f%%", pct_overlap), "")),
      hjust  = -0.15,
      size   = 3.2,
      colour = "grey30"
    ) +
    scale_fill_manual(values = FACTOR_COLOURS[ds_label], name = NULL) +
    scale_x_continuous(limits = c(0, 125),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "% overlap", y = NULL) +
    base_theme +
    theme(
      axis.text.y     = if (show_y) element_text(size = 12) else element_blank(),
      axis.ticks.y    = if (show_y) element_line() else element_blank(),
      legend.position = "none"
    )
}

# Assemble one row (one direction) across all 4 datasets ---------------------
# Column order per dataset: counts, pct
# Width ratio per pair: 3 (counts) : 1 (pct)
# y-axis labels only on the very first column

make_row <- function(dir_label) {
  col <- ifelse(dir_label == "Higher in AD", "#C0392B", "#2980B9")
  pairs <- lapply(seq_along(DATASETS), function(i) {
    ds       <- DATASETS[i]
    # celltype labels only on the first counts panel; direction label sits left of it
    show_y   <- (i == 1)
    show_ttl <- TRUE

    counts_p <- make_counts_col(ds, dir_label,
                                show_y     = show_y,
                                show_title = show_ttl)
    pct_p    <- make_pct_col(ds, dir_label, show_y = FALSE)

    counts_p + pct_p + plot_layout(widths = c(3, 1))
  })

  # Rotated direction label — sits immediately left of the celltype labels
  label_panel <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = dir_label,
             angle = 90, size = 5, fontface = "bold", colour = col,
             hjust = 0.5, vjust = 0.5) +
    theme_void()

  content <- patchwork::wrap_plots(pairs, nrow = 1)

  # narrow direction label | celltype-labelled content
  (label_panel | content) + plot_layout(widths = c(0.04, 1))
}

row_higher <- make_row("Higher in AD")
row_lower  <- make_row("Lower in AD")

# Shared legend
shared_legend <- cowplot::get_legend(
  ggplot(plot_data, aes(x = n_genes, y = celltype, fill = factor)) +
    geom_col() +
    fill_scale +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.5, "cm"),
          legend.text     = element_text(size = 12))
)

# Final assembly: rows stacked, no label column
fig1 <- (row_higher / row_lower / patchwork::wrap_elements(shared_legend)) +
  plot_layout(heights = c(1, 1, 0.08)) +
  plot_annotation(
    title    = "Gene counts and dCoPA overlap by celltype",
    subtitle = "Each dataset pair: left = # genes (DE dataset vs dcopa), right = % of dcopa genes recovered",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 12, colour = "grey40")
    )
  )

ggsave(file.path(SAVE_DIR, "fig1_counts_and_pct_combined.pdf"),
       fig1, width = FIG1_W, height = FIG1_H)
message("  Saved: fig1_counts_and_pct_combined.pdf")


# ============================================================
# FIGURE 2 - MIRROR BAR: higher vs lower, all 5 factors dodged
# ============================================================
message("Building Figure 2: mirror bar ...")

mirror_data <- plot_data |>
  dplyr::mutate(
    n_plot = ifelse(direction == "Higher in AD", n_genes, -n_genes)
  )

x_max <- max(abs(mirror_data$n_plot), na.rm = TRUE) * 1.12

fig2 <- ggplot(mirror_data,
               aes(x = n_plot, y = celltype, fill = factor)) +
  geom_col(position = position_dodge(width = 0.85),
           width = 0.8, colour = NA) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey20") +
  fill_scale +
  scale_x_continuous(
    limits = c(-x_max, x_max),
    labels = function(x) abs(x),
    expand = expansion(mult = 0.02)
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(
    title    = "Gene counts: higher vs lower in AD",
    subtitle = "Right = higher in AD; Left = lower in AD",
    x        = "Number of genes",
    y        = NULL
  ) +
  base_theme +
  theme(plot.margin = margin(1, 2, 1, 1, "cm"))

ggsave(file.path(SAVE_DIR, "fig2_mirror_bar.pdf"),
       fig2, width = FIG2_W, height = FIG2_H)
message("  Saved: fig2_mirror_bar.pdf")


message("\n====== All figures saved to ", SAVE_DIR, " ======")
message("  fig1 - 8-column layout: counts + % overlap for each DE dataset vs dcopa, by direction")
message("  fig2 - Mirror bar: higher (right) vs lower (left) for all 5 factors")

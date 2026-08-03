# Fig. S11 — standalone reproduction of the dCoPA-overlap significance heatmaps.
#
# Cross-celltype Fisher's exact tests for dCoPA gene-set overlap (AD modules vs
# CTRL modules), rendered as celltype x celltype -log10(FDR) heatmaps for the
# "lower in more severe" direction, for both MTG and DFC, plus a shared legend.
#
# This is a self-contained copy of fig_7/v7.1/panel_C_D_fisher_dcopa_combined.R.
# It REFERENCES the precomputed dCoPA gene lists and the bulk-gene universe that
# already live under fig_7/ (raw data is NOT regenerated here); only the Fisher
# tests and the figures are recomputed. Outputs are written next to this script.
#
# Outputs (in fig_s11/v3, all vector SVG — assemble_figure.py composes the figure):
#   panel_a_MTG.svg        (left panel; Fig. S11a)
#   panel_b_DFC.svg        (right panel; Fig. S11b)
#   legend.svg             (shared -log10(FDR) legend, transparent background)
#   panel_{a_MTG,b_DFC}_padj.csv   (BH-adjusted p-value matrices)

library(data.table)
library(tidyverse)
library(showtext)
library(cowplot)
library(svglite)
library(ggpubr)

# Fonts: render everything in Arial. Liberation Sans is metric-compatible with
# Arial and ships on this system; register it under the name "Arial".
.LIBSANS <- "/usr/share/fonts/truetype/liberation"
if (file.exists(file.path(.LIBSANS, "LiberationSans-Regular.ttf"))) {
  sysfonts::font_add("Arial",
    regular    = file.path(.LIBSANS, "LiberationSans-Regular.ttf"),
    bold       = file.path(.LIBSANS, "LiberationSans-Bold.ttf"),
    italic     = file.path(.LIBSANS, "LiberationSans-Italic.ttf"),
    bolditalic = file.path(.LIBSANS, "LiberationSans-BoldItalic.ttf"))
  FONT <- "Arial"
} else {
  FONT <- "sans"   # fallback: device default sans-serif
}
showtext_auto()

# ============================================================
# USER INPUT
# ============================================================

# Shared settings (can be overridden per run via fdr_thresh / neglog_cap fields)
MIN_GENES   <- 5
FDR_THRESH  <- 0.05
NEGLOG_CAP  <- 7

# Only the "lower in more severe" direction is featured in Fig. S11 (matches the
# panel_lower reference). Set to c("higher", "lower") to also emit the higher panel.
DIRECTIONS  <- c("lower")

# Fill scale: two-colour blue -> red ramp over the full -log10(FDR) range
# (low significance reads blue, high significance reads red).
FILL_SCALE <- ggplot2::scale_fill_gradient(
  low = "#2166AC", high = "#B2182B",
  limits = c(0, NEGLOG_CAP), name = "-log10(FDR)")

# Raw data lives in fig_7/ — reference it, do not regenerate.
FIG7        <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7")
UNIVERSE    <- file.path(FIG7, "v3", "all_bulk_genes.txt")
OUT_DIR     <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s11/v3")

CT_RENAME <- c(
  "Lamp5"     = "LAMP5",
  "Lamp5 Lhx6"= "LAMP5 LHX6",
  "Pax6"      = "PAX6",
  "Pvalb"     = "PVALB",
  "Sncg"      = "SNCG",
  "Sst"       = "SST",
  "Vip"       = "VIP",
  "Sst Chodl" = "SST CHODL"
)

CT_ORDER_REQUESTED <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP",
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "SST CHODL", "SST", "Chandelier",
              "PVALB", "LAMP5", "LAMP5 LHX6", "PAX6",  "SNCG",
              "VIP",  "Astrocyte", "Oligodendrocyte", "OPC", "Microglia-PVM", "Endothelial", "VLMC")

# Runs — set_keys are read automatically from the first column of dcopa_path.
RUNS <- list(
  list(
    dcopa_path    = file.path(FIG7, "v7.1", "panel_C_D_MTG_dcopa_genelist.csv"),
    save_dir      = OUT_DIR,
    universe_path = UNIVERSE,
    file_prefix   = "panel_a_MTG"
  ),
  list(
    dcopa_path    = file.path(FIG7, "v7.1", "panel_C_D_DFC_dcopa_genelist.csv"),
    save_dir      = OUT_DIR,
    universe_path = UNIVERSE,
    file_prefix   = "panel_b_DFC"
  )
)

# ============================================================
# HELPERS
# ============================================================
# Save individual panels as vector SVG only (the assembly composes a vector figure
# from these; a duplicate PDF per panel is not needed).
save_plot <- function(path_no_ext, plot_expr, width, height, svg_bg = "white") {
  svglite(paste0(path_no_ext, ".svg"), width = width, height = height, bg = svg_bg)
  print(plot_expr)
  dev.off()
  message("Saved: ", basename(path_no_ext), ".svg")
}

# ============================================================
# ANALYSIS FUNCTION
# Returns: full_table + set_keys
# ============================================================
run_fisher_analysis <- function(dcopa_path,
                                save_dir,
                                universe_path = NULL,
                                min_genes     = MIN_GENES,
                                file_prefix   = "panel_C_D") {

  message("\n=== Fisher analysis: ", basename(dcopa_path), " ===")
  dcopa_raw <- fread(dcopa_path, data.table = FALSE)

  # Derive set_keys from unique values of the first column
  set_keys <- unique(dcopa_raw[[1]])
  message("  set_keys detected: ", paste(set_keys, collapse = " | "))

  dcopa <- dcopa_raw |>
    dplyr::mutate(
      Celltype  = dplyr::recode(Celltype, !!!CT_RENAME),
      direction = dplyr::case_when(
        Direction == "Higher in more severe" ~ "higher",
        Direction == "Lower in more severe"  ~ "lower",
        TRUE                                 ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(direction), Comparison %in% set_keys)

  message("  ", nrow(dcopa), " rows retained after normalisation")

  all_cts  <- sort(unique(dcopa$Celltype))
  all_dirs <- c("higher", "lower")

  gene_sets <- setNames(lapply(set_keys, function(comp) {
    sub <- dcopa[dcopa$Comparison == comp, ]
    setNames(lapply(all_cts, function(ct) {
      setNames(lapply(all_dirs, function(dir) {
        unique(sub$Gene[sub$Celltype == ct & sub$direction == dir])
      }), all_dirs)
    }), all_cts)
  }), set_keys)

  if (!is.null(universe_path)) {
    user_universe  <- unique(readLines(universe_path))
    user_universe  <- user_universe[nzchar(user_universe)]
    message("  User-supplied universe: ", length(user_universe), " genes")
    global_universe <- setNames(lapply(all_dirs, function(d) user_universe), all_dirs)
  } else {
    message("  No universe_path — using union of all groups × all celltypes per direction")
    global_universe <- setNames(lapply(all_dirs, function(dir) {
      unique(unlist(lapply(set_keys, function(comp) {
        unlist(lapply(all_cts, function(ct) gene_sets[[comp]][[ct]][[dir]]))
      })))
    }), all_dirs)
  }

  run_fisher <- function(genes_A, genes_B, universe) {
    n_univ <- length(universe)
    if (n_univ == 0) return(NULL)
    if (length(genes_A) < min_genes || length(genes_B) < min_genes) return(NULL)
    a <- length(intersect(genes_A, genes_B))
    b <- length(genes_A) - a
    c <- length(genes_B) - a
    d <- n_univ - a - b - c
    if (d < 0) { warning("Universe smaller than gene sets; skipping."); return(NULL) }
    ft <- fisher.test(matrix(c(a, c, b, d), nrow = 2), alternative = "greater")
    data.frame(n_G1 = length(genes_A), n_G2 = length(genes_B),
               n_overlap = a, n_universe = n_univ,
               odds_ratio = unname(ft$estimate), p_value = ft$p.value,
               stringsAsFactors = FALSE)
  }

  pairs <- combn(set_keys, 2, simplify = FALSE)
  message("Running cross-celltype Fisher's tests (",
          length(pairs), " pair(s) × ", length(all_cts), "² celltypes × ",
          length(all_dirs), " directions) ...")

  results <- lapply(pairs, function(pair) {
    g1 <- pair[1]; g2 <- pair[2]
    rows <- lapply(all_cts, function(ct1) {
      lapply(all_cts, function(ct2) {
        lapply(all_dirs, function(dir) {
          res <- run_fisher(gene_sets[[g1]][[ct1]][[dir]],
                            gene_sets[[g2]][[ct2]][[dir]],
                            global_universe[[dir]])
          if (is.null(res)) return(NULL)
          cbind(data.frame(Group1 = g1, Group2 = g2,
                           Celltype1 = ct1, Celltype2 = ct2,
                           Direction = dir, stringsAsFactors = FALSE), res)
        })
      })
    })
    dplyr::bind_rows(Filter(Negate(is.null), unlist(rows, recursive = FALSE)))
  })

  full_table <- dplyr::bind_rows(results)
  full_table$p_adj <- p.adjust(full_table$p_value, method = "BH")
  message("  ", nrow(full_table), " tests performed, ",
          sum(full_table$p_adj < 0.05, na.rm = TRUE), " significant at FDR < 0.05")

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  list(full_table = full_table, set_keys = set_keys)
}

# ============================================================
# VISUALISATION FUNCTION
# ============================================================
run_fisher_viz <- function(full_table,
                           save_dir,
                           set_keys,
                           axis_label         = NULL,
                           ct_order_requested = CT_ORDER_REQUESTED,
                           fdr_thresh         = FDR_THRESH,
                           neglog_cap         = NEGLOG_CAP,
                           directions         = DIRECTIONS,
                           file_prefix        = "panel_C_D") {

  message("\n=== Fisher visualisation ===")

  if (is.null(axis_label)) axis_label <- setNames(set_keys, set_keys)

  ft <- full_table |>
    dplyr::mutate(
      neglog10_padj = pmin(-log10(p_adj + 1e-300), neglog_cap),
      sig           = p_adj < fdr_thresh
    )

  data_cts <- union(unique(ft$Celltype1), unique(ft$Celltype2))
  ct_order <- c(ct_order_requested,
                setdiff(data_cts, ct_order_requested))
  n_cts <- length(ct_order)

  plot_df <- ft |>
    dplyr::mutate(ct1 = factor(Celltype1, levels = ct_order),
                  ct2 = factor(Celltype2, levels = ct_order)) |>
    tidyr::complete(ct1, ct2, Direction) |>
    dplyr::mutate(is_tested = !is.na(neglog10_padj),
                  sig       = dplyr::coalesce(sig, FALSE))

  sig_df <- dplyr::filter(plot_df, sig)
  na_df  <- dplyr::filter(plot_df, !is_tested)

  build_ct_heatmap <- function(df, sig_pts, na_pts,
                               group1_label, group2_label,
                               show_legend = TRUE) {
    ggplot(df, aes(x = ct2, y = ct1)) +
      geom_tile(data = na_pts, fill = "white", colour = "grey80", linewidth = 0.3) +
      geom_text(data = na_pts, aes(label = "n/a"), family = FONT,
                colour = "grey70", size = 6, inherit.aes = TRUE) +
      geom_tile(data = dplyr::filter(df, is_tested),
                aes(fill = neglog10_padj), colour = "white", linewidth = 0.3) +
      geom_text(data = sig_pts, aes(label = "*"), family = FONT,
                colour = "white", size = 14, vjust = 0.75, inherit.aes = TRUE) +
      # Mark the diagonal (same celltype on both axes) with a dashed black line.
      geom_segment(data = data.frame(x = 0.5, y = n_cts + 0.5,
                                     xend = n_cts + 0.5, yend = 0.5),
                   aes(x = x, y = y, xend = xend, yend = yend),
                   inherit.aes = FALSE, colour = "black",
                   linewidth = 0.9, linetype = "22") +
      FILL_SCALE +
      coord_fixed() +
      scale_x_discrete(drop = FALSE, position = "bottom") +
      scale_y_discrete(drop = FALSE, limits = rev(ct_order)) +
      theme_bw(base_size = 26, base_family = FONT) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 24),
        axis.text.y     = element_text(size = 24),
        axis.title.x    = element_text(size = 26, face = "plain", margin = margin(t = 6)),
        axis.title.y    = element_text(size = 26, face = "plain", margin = margin(r = 6)),
        legend.title    = element_text(size = 24),
        legend.text     = element_text(size = 22),
        legend.position = if (show_legend) "right" else "none",
        plot.title      = element_blank(),
        plot.subtitle   = element_blank()
      ) +
      labs(x = group2_label, y = group1_label)
  }

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  g1_label <- axis_label[set_keys[1]]
  g2_label <- axis_label[set_keys[2]]
  plot_w   <- n_cts * 0.5 + 4
  plot_h   <- n_cts * 0.5 + 3

  p_for_legend <- NULL

  for (dir_name in directions) {
    panel_name <- dir_name
    df_dir  <- dplyr::filter(plot_df, Direction == dir_name)
    sig_dir <- dplyr::filter(sig_df,  Direction == dir_name)
    na_dir  <- dplyr::filter(na_df,   Direction == dir_name)

    if (nrow(df_dir) == 0) {
      message("Skipping '", dir_name, "' direction — no data"); next
    }

    p_noleg <- build_ct_heatmap(df_dir, sig_dir, na_dir, g1_label, g2_label,
                                show_legend = FALSE)
    # file_prefix already reflects the panel (panel_a_MTG / panel_b_DFC)
    save_plot(file.path(save_dir, file_prefix), p_noleg, plot_w, plot_h)

    # Save the adjusted p-values as a matrix, laid out as in the heatmap
    padj_wide <- df_dir |>
      dplyr::select(ct1, ct2, p_adj) |>
      tidyr::pivot_wider(names_from = ct2, values_from = p_adj, names_sort = TRUE) |>
      dplyr::arrange(ct1)
    csv_path <- file.path(save_dir, paste0(file_prefix, "_padj.csv"))
    data.table::fwrite(padj_wide, csv_path)
    message("Saved: ", basename(csv_path))

    if (is.null(p_for_legend))
      p_for_legend <- build_ct_heatmap(df_dir, sig_dir, na_dir, g1_label, g2_label,
                                       show_legend = TRUE)
  }
  if (!is.null(p_for_legend)) {
    legend_grob  <- cowplot::get_legend(
      p_for_legend + theme(
        legend.text          = element_text(size = 16),
        legend.title         = element_text(size = 18),
        legend.justification = "center",
        legend.title.align   = 0.5
      )
    )
    # Transparent background so the legend never occludes a panel when overlaid.
    p_legend_out <- ggpubr::as_ggplot(legend_grob) +
      theme(plot.margin     = margin(0, 0, 0, 0),
            plot.background  = element_rect(fill = NA, colour = NA),
            panel.background = element_rect(fill = NA, colour = NA))
    # One shared legend for both panels (the -log10(FDR) scale is identical).
    save_plot(file.path(save_dir, "legend"),
              p_legend_out, 4, 3, svg_bg = "transparent")
  }

  invisible(NULL)
}

# ============================================================
# MAIN LOOP
# ============================================================
for (run in RUNS) {
  result <- run_fisher_analysis(
    dcopa_path    = run$dcopa_path,
    save_dir      = run$save_dir,
    universe_path = run$universe_path,
    file_prefix   = run$file_prefix %||% "panel_C_D"
  )

  run_fisher_viz(
    full_table  = result$full_table,
    save_dir    = run$save_dir,
    set_keys    = result$set_keys,
    axis_label  = run$axis_label,   # NULL → use set_keys as-is
    file_prefix = run$file_prefix %||% "panel_C_D"
  )
}

message("\nAll runs complete.")

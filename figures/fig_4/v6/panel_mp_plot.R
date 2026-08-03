# panel_mp_plot.R
# Shared builders for Fig. 4 panels m-p (per-module projection barplots) and the
# nested dataset-correlation brackets overlaid on them. Sourced by both fig_4_v6.R
# (full figure) and the standalone panel_mp.R, so this styling lives in one place.
# Requires (loaded by the caller): ggplot2, cowplot, dplyr, scales.

DATASET_LABELS <- c(
  "Jorstad et al. 2023 (normal human dorsolateral prefrontal cortex)",
  "Gabitto et al. 2024 (normal human dorsolateral prefrontal cortex)",
  "Liu et al. 2025 (normal human dorsolateral prefrontal cortex)",
  "Morabito et al. 2021 (normal human dorsolateral prefrontal cortex)")

# Build one projection barplot for module i. `dfs` / `ses` are length-4 lists of
# module x celltype mean / SE data frames in dataset order (Jorstad, Gabitto, Liu,
# Morabito). Returns list(plot = ggplot, cormat = 4x4 pairwise-cor matrix).
make_proj_panel <- function(i, dfs, ses, allcts_cap, proj_pal, xlab) {
  n <- length(allcts_cap)
  plotdf <- data.frame(
    ct      = rep(allcts_cap, 4),
    dataset = rep(DATASET_LABELS, each = n),
    ind     = c(unlist(dfs[[1]][i, ]), unlist(dfs[[2]][i, ]), unlist(dfs[[3]][i, ]), unlist(dfs[[4]][i, ])),
    ind_se  = c(unlist(ses[[1]][i, ]), unlist(ses[[2]][i, ]), unlist(ses[[3]][i, ]), unlist(ses[[4]][i, ]))
  ) |>
    dplyr::mutate(dataset = factor(dataset, levels = DATASET_LABELS),
                  ct      = factor(ct, levels = allcts_cap))

  tempdf <- data.frame("Jorstad 2023" = unlist(dfs[[1]][i, ]),
                       "Gabitto 2024" = unlist(dfs[[2]][i, ]),
                       "Liu 2025"     = unlist(dfs[[3]][i, ]),
                       "Morabito 2021"= unlist(dfs[[4]][i, ]), check.names = FALSE)
  cormat <- cor(tempdf, use = "pairwise.complete.obs")

  p <- ggplot(plotdf, aes(x = ct, y = ind, fill = ct)) +
    theme_classic() +
    geom_col(position = position_dodge(), alpha = 0.8) +
    geom_errorbar(aes(ymin = ind - 2 * ind_se, ymax = ind + 2 * ind_se),
                  width = 0.2, linewidth = 0.3, position = position_dodge(0.5)) +
    scale_fill_manual(values = proj_pal, na.value = "grey70") +
    theme(text = element_text(family = "sans", color = "black", size = 12),
          legend.position = "none",
          axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
          axis.text.x  = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
          strip.text   = element_text(color = "black"),
          strip.background = element_rect(fill = "white")) +
    facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free_y") +
    labs(y = xlab, x = "") +
    scale_y_continuous(breaks = c(0, 0.5, 1))

  list(plot = p, cormat = cormat)
}

# Overlay nested dataset-correlation brackets + "r = 0.XX" labels, in the style of
# fig_3/v4's add_corr_brackets(). One bracket per dataset pair (C(4,2)=6), nested.
# The bracket geometry was recovered from the hand-drawn Inkscape annotation
# panel_1_4_brackets.svg (layer g1-9): each entry's `pts` are the bracket's 4 corner
# points (serif-top, spine-top, spine-bottom, serif-bottom) as ABSOLUTE positions in
# the reference's mm frame, where the base 5.5 x 5 in (396 x 360 pt) panel sits at the
# origin and the brackets occupy the right margin. The panel keeps its base width; the
# extra (width_in - 5.5) is blank margin holding the annotation. `cormat` is the 4x4
# correlation matrix in dataset order (Jorstad, Gabitto, Liu, Morabito).
# Returns a cowplot ggdraw object; save it with width = width_in, height = height_in.
add_corr_brackets <- function(p, cormat, width_in = 6.6, height_in = 5) {
  mm2pt <- 72 / 25.4
  W <- width_in * 72; H <- height_in * 72
  right_margin_pt <- W - 5.5 * 72 + 5.5     # keep the 5.5in panel, blank to the right
  gap <- 3; font <- 10                      # label gap (pt) and size
  bl <- list(  # pts (mm): serif-top, spine-top, spine-bottom, serif-bottom; i,j = pair
    list(pts = list(c(139.2, 15.9), c(140.6, 15.9), c(140.6, 38.6), c(139.1, 38.6)), i = 1, j = 2),
    list(pts = list(c(139.1, 40.7), c(140.6, 40.7), c(140.6, 62.7), c(139.0, 62.7)), i = 2, j = 3),
    list(pts = list(c(139.0, 64.8), c(140.4, 64.8), c(140.4, 85.6), c(138.9, 85.6)), i = 3, j = 4),
    list(pts = list(c(144.9, 15.9), c(146.2, 15.9), c(146.2, 62.7), c(144.8, 62.7)), i = 1, j = 3),
    list(pts = list(c(150.9, 39.6), c(152.2, 39.6), c(152.2, 85.6), c(150.8, 85.6)), i = 2, j = 4),
    list(pts = list(c(156.2, 15.9), c(157.5, 15.9), c(157.5, 85.7), c(156.1, 85.7)), i = 1, j = 4))
  xn <- function(xpt) xpt / W
  yn <- function(ypt) 1 - ypt / H
  g <- ggdraw(p + theme(plot.margin = margin(5.5, right_margin_pt, 5.5, 5.5, "pt")))
  for (b in bl) {
    px <- vapply(b$pts, function(q) q[1] * mm2pt, numeric(1))
    py <- vapply(b$pts, function(q) q[2] * mm2pt, numeric(1))
    g <- g + draw_line(xn(px), yn(py), color = "black", linewidth = 0.4)
    spine_x    <- b$pts[[2]][1] * mm2pt
    spine_ymid <- (b$pts[[2]][2] + b$pts[[3]][2]) / 2 * mm2pt
    # truncate toward zero so a near-perfect r (e.g. 0.996) shows "0.99", not "1.00"
    r_disp <- trunc(cormat[b$i, b$j] * 100 + sign(cormat[b$i, b$j]) * 1e-9) / 100
    g <- g + draw_label(sprintf("r = %.2f", r_disp),
                        x = xn(spine_x + gap + font / 2), y = yn(spine_ymid),
                        angle = 270, size = font, fontfamily = "sans")
  }
  g
}

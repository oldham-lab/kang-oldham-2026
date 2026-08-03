# panel_S_plot.R
# Shared renderer for Fig. 4 panel S (projection-correlation heatmap).
# Sourced by both fig_4_v6.R (full figure) and panel_S.R (fast standalone), so the
# panel-S styling lives in exactly one place.
# Requires (loaded by the caller): ggplot2, cowplot, ggdendro, scales.

make_panel_S <- function(cormeans, save_path) {
  mat <- cormeans
  row_dend <- as.dendrogram(hclust(dist(mat)))
  col_dend <- as.dendrogram(hclust(dist(t(mat))))
  mat_ord  <- mat[order.dendrogram(row_dend), order.dendrogram(col_dend)]

  cor_df <- as.data.frame(as.table(mat_ord))
  colnames(cor_df) <- c("row_lab", "col_lab", "value")
  cor_df$row_lab <- factor(cor_df$row_lab, levels = rownames(mat_ord))
  cor_df$col_lab <- factor(cor_df$col_lab, levels = colnames(mat_ord))

  hmap <- ggplot(cor_df, aes(x = col_lab, y = row_lab, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(value, 2)), size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1),
                         name = "Correlation") +
    coord_fixed() +
    theme_minimal() +
    theme(axis.text.x           = element_blank(),
          axis.ticks.x          = element_blank(),
          axis.text.y           = element_text(size = 13, hjust = 1),
          axis.title            = element_blank(),
          legend.position       = "bottom",
          legend.direction      = "horizontal",
          legend.title.position = "top",
          legend.key.height     = unit(0.4, "cm"),
          legend.key.width      = unit(0.6, "cm"),
          legend.title          = element_text(size = 13, hjust = 0.5),
          legend.text           = element_text(size = 13),
          plot.margin           = margin(2, 2, 2, 2))

  col_segs   <- segment(dendro_data(col_dend, type = "rectangle"))
  col_dend_p <- ggplot(col_segs) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.4) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_void()

  row_segs   <- segment(dendro_data(row_dend, type = "rectangle"))
  row_dend_p <- ggplot(row_segs) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.4) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    theme_void()

  # Heatmap block (dendrograms on top/right, no legend); the legend and the overall-
  # mean annotation are stacked in a right-hand column so the legend sits beneath the
  # "Overall mean" text rather than under the heatmap.
  hmap_nl <- hmap + theme(legend.position = "none")
  p1 <- insert_xaxis_grob(hmap_nl, col_dend_p, grid::unit(0.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1,      row_dend_p, grid::unit(0.2, "null"), position = "right")
  leg <- get_legend(hmap)
  ov  <- mean(cormeans[upper.tri(cormeans)])
  p <- ggdraw() +
    draw_plot(p2, x = 0, y = 0, width = 0.63, height = 1) +
    draw_label(sprintf("Overall mean:\nr = %.2f", ov),
               x = 0.82, y = 0.62, size = 13, hjust = 0.5, lineheight = 0.95) +
    draw_grob(leg, x = 0.63, y = 0.16, width = 0.37, height = 0.30)

  ggsave(p, file = save_path, height = 3, width = 5.2)
  invisible(p)
}

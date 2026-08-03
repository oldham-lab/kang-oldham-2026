# ============================================================
# Approximate Dendrogram33 using ComplexHeatmap
# Dendrogram33 referring to the metamodule dendrogram (consensusMin, Jorstad) that is labeled
# Required packages: ComplexHeatmap, circlize, grid
#   install via: BiocManager::install("ComplexHeatmap"); install.packages("circlize")
#
# Required input files (update paths as needed):
#   - mod_eig.csv                            : 24 x 101 eigengene matrix
#   - branchpoint_table_modeig_with_genes.csv: branchpoint annotations
#   - indiv_leaf_modgene_per.csv             : per-module mod_per and gene_per
# Output:
#   - dendrogram33.pdf
# ============================================================

library(ComplexHeatmap)
library(circlize)
library(grid)

# ---- File paths ----
MOD_EIG_FILE  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/mod_eig.csv")
BP_GENES_FILE <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/branchpoint_table_modeig_with_genes.csv")
LEAF_PER_FILE <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/indiv_leaf_modgene_per.csv")
OUTPUT_FILE   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/dendrogram33_R.pdf")

# ---- 1. Load data ----
eig      <- read.csv(MOD_EIG_FILE,  header = TRUE, row.names = 1, check.names = FALSE)
modules  <- colnames(eig)
n        <- length(modules)
eig_mat  <- as.matrix(eig)

bp_df    <- read.csv(BP_GENES_FILE, header = TRUE, check.names = FALSE)
leaf_per <- read.csv(LEAF_PER_FILE, header = TRUE, row.names = 1, check.names = FALSE)

# ---- 2. Pearson correlation + complete-linkage clustering ----
corr_mat  <- cor(eig_mat, method = "pearson")
dist_mat  <- as.dist(1 - corr_mat)
hc        <- hclust(dist_mat, method = "complete")
cut_height <- 0.3

# ---- 3. Identify branchpoints above cut ----
node_heights <- hc$height
merge        <- hc$merge

# Iterative x-position computation (avoids stack overflow)
leaf_pos   <- setNames(seq_along(hc$order) - 0.5, hc$order)
node_x_pos <- numeric(n - 1)
for (i in seq_len(n - 1)) {
  left  <- merge[i, 1]
  right <- merge[i, 2]
  xl <- if (left  > 0) leaf_pos[left]  else node_x_pos[-left]
  xr <- if (right > 0) leaf_pos[right] else node_x_pos[-right]
  node_x_pos[i] <- (xl + xr) / 2
}

# Branchpoints sorted by height (matching Python BP labeling)
bp_indices <- which(node_heights >= cut_height)
bp_indices <- bp_indices[order(node_heights[bp_indices])]

# Build branchpoint annotation list
pct_mods_lookup  <- setNames(bp_df[["Pct_of_Total"]],       bp_df[["Branchpoint_ID"]])
pct_genes_lookup <- setNames(bp_df[["Pct_of_Total_genes"]], bp_df[["Branchpoint_ID"]])

bp_list <- lapply(seq_along(bp_indices), function(j) {
  i     <- bp_indices[j]
  bp_id <- sprintf("BP%02d", j)
  list(
    merge_row  = i,
    height     = node_heights[i],
    bp_id      = bp_id,
    pct_mods   = pct_mods_lookup[[bp_id]],
    pct_genes  = pct_genes_lookup[[bp_id]]
  )
})

# ---- 4. Per-module bar annotations ----
# Reorder leaf_per to match hclust leaf order
ordered_modules <- modules[hc$order]
mod_per_vals    <- leaf_per[ordered_modules, "mod_per"]
gene_per_vals   <- leaf_per[ordered_modules, "gene_per"]

# Bar annotation function
make_bar_anno <- function(values, color, label) {
  anno_df <- data.frame(x = values)
  HeatmapAnnotation(
    df  = anno_df,
    col = list(x = colorRamp2(c(0, max(values)), c("white", color))),
    annotation_name_side = "left",
    annotation_label     = label,
    annotation_name_gp   = gpar(fontsize = 8),
    show_legend          = FALSE,
    bar                  = anno_barplot(
      values,
      gp         = gpar(fill = color, col = NA),
      axis       = FALSE,
      bar_width  = 1,
      add_numbers = TRUE,
      numbers_gp  = gpar(fontsize = 7)
    ),
    which = "column",
    height = unit(1.2, "cm")
  )
}

mod_bar_anno  <- columnAnnotation(
  `Mod %` = anno_barplot(
    mod_per_vals,
    gp          = gpar(fill = "#4878cf", col = NA),
    axis        = FALSE,
    bar_width   = 1,
    add_numbers = TRUE,
    numbers_gp  = gpar(fontsize = 7)
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 8),
  height = unit(1.2, "cm")
)

gene_bar_anno <- columnAnnotation(
  `Gene %` = anno_barplot(
    gene_per_vals,
    gp          = gpar(fill = "#6acc65", col = NA),
    axis        = FALSE,
    bar_width   = 1,
    add_numbers = TRUE,
    numbers_gp  = gpar(fontsize = 7)
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 8),
  height = unit(1.2, "cm")
)

combined_anno <- c(mod_bar_anno, gene_bar_anno)

# ---- 5. Heatmap color scale ----
eig_ordered <- eig_mat[, ordered_modules]
vmax        <- quantile(abs(eig_ordered), 0.98)
col_fun     <- colorRamp2(c(-vmax, 0, vmax), c("blue", "white", "red"))

# ---- 6. Draw heatmap ----
ht <- Heatmap(
  eig_ordered,
  name                   = "Eigengene",
  col                    = col_fun,
  cluster_rows           = FALSE,
  cluster_columns        = as.dendrogram(hc),
  column_dend_height     = unit(6, "cm"),
  column_dend_gp         = gpar(lwd = 3),
  show_column_names      = TRUE,
  column_names_rot       = 90,
  column_names_gp        = gpar(fontsize = 12),
  row_names_gp           = gpar(fontsize = 13),
  row_title              = "Cell type",
  row_title_gp           = gpar(fontsize = 16),
  column_title           = "Jorstad metamodules, consensusMin",
  column_title_gp        = gpar(fontsize = 15),
  column_title_side      = "top",
  top_annotation         = combined_anno,
  heatmap_legend_param   = list(
    title          = "Eigengene value",
    direction      = "horizontal",
    title_position = "topcenter",
    legend_width   = unit(4, "cm"),
    labels_gp      = gpar(fontsize = 9),
    title_gp       = gpar(fontsize = 11)
  )
)

# ---- 7. Draw and annotate branchpoints ----
pdf(OUTPUT_FILE, width = 40, height = 24)

draw(ht,
     heatmap_legend_side = "bottom",
     padding             = unit(c(15, 2, 2, 2), "mm"))

decorate_column_dend("Eigengene", {
  max_h   <- max(hc$height)
  x_scale <- 1 / n
  y_scale <- 1 / max_h

  # Cut line
  grid.lines(
    x  = unit(c(0, 1), "npc"),
    y  = unit(rep(cut_height / max_h, 2), "npc"),
    gp = gpar(col = "red", lty = 2, lwd = 2.5)
  )

  for (bp in bp_list) {
    x_npc <- node_x_pos[bp$merge_row] * x_scale
    y_npc <- bp$height * y_scale

    label_text <- sprintf("%s\n%.1f%% (m)\n%.1f%% (g)",
                          bp$bp_id, bp$pct_mods, bp$pct_genes)

    # Background box
    grid.rect(
      x     = unit(x_npc, "npc"),
      y     = unit(y_npc, "npc"),
      width = unit(1.8, "cm"), height = unit(1.1, "cm"),
      just  = c("center", "bottom"),
      gp    = gpar(fill = rgb(1, 1, 0.9, 0.5), col = NA)
    )
    # Label text
    grid.text(
      label = label_text,
      x     = unit(x_npc, "npc"),
      y     = unit(y_npc, "npc"),
      just  = c("center", "bottom"),
      gp    = gpar(fontsize = 9, fontface = "bold",
                   col = rgb(0.55, 0, 0, 0.75))
    )
  }
})

dev.off()
message("Saved: ", OUTPUT_FILE)

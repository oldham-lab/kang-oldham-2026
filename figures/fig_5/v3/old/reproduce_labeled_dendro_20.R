# ============================================================
# Reproduce Dendrogram20 using ComplexHeatmap
# Inputs required (update paths as needed):
#   - test.csv                          : 96x96 correlation matrix
#   - mod_eig1.csv                      : 24x96 eigengene matrix
#   - branchpoint_table3_with_genes.csv : branchpoint table with % mods and % genes
# ============================================================

# Produced by Claude
# The branchpoint labels produced by this code are all over the place, needs to be fixed

library(ComplexHeatmap)
library(circlize)
library(grid)

# ---- 1. Load data ----
corr_df   <- read.csv("/home/gugene/test/test.csv", header = TRUE, check.names = FALSE)
eig_df    <- read.csv(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/mod_eig1.csv"), header = TRUE, check.names = FALSE)
bp_df     <- read.csv(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/branchpoint_table3_with_genes.csv"), header = TRUE,
                      check.names = FALSE)


modules   <- colnames(corr_df)          # 96 module names
corr_mat  <- as.matrix(corr_df)
rownames(corr_mat) <- modules

eig_mat   <- as.matrix(eig_df)          # 24 x 96
rownames(eig_mat) <- paste0("S", seq_len(nrow(eig_mat)))

# ---- 2. Complete-linkage clustering on distance = 1 - r ----
dist_mat  <- as.dist(1 - corr_mat)
hc        <- hclust(dist_mat, method = "complete")
cut_height <- 0.3

# ---- 3. Build branchpoint annotation labels ----
# For every node above the cut, compose a 3-line label:
#   BP##  /  X.X% (m)  /  X.X% (g)
# We'll draw these as text grobs on top of the dendrogram.

# Helper: get all leaf indices under a node in hclust merge matrix
get_leaves <- function(node, merge, n) {
  if (node <= n) return(node)   # positive index = leaf
  row <- -node                  # negative index in merge = internal node
  c(get_leaves(merge[row, 1], merge, n),
    get_leaves(merge[row, 2], merge, n))
}

n      <- length(modules)
merge  <- hc$merge          # (n-1) x 2 matrix

# Identify internal nodes above cut_height
node_heights <- hc$height   # length n-1

# Branchpoints: internal nodes with height >= cut_height
bp_indices <- which(node_heights >= cut_height)  # row indices in merge/height

# Sort by height ascending (matching the Python BP labeling)
bp_indices <- bp_indices[order(node_heights[bp_indices])]

# Build lookup from BP ID to pct values
bp_lookup <- setNames(
  list(),
  character(0)
)
for (i in seq_len(nrow(bp_df))) {
  bp_lookup[[ bp_df[i, "Branchpoint_ID"] ]] <- list(
    pct_mods  = bp_df[i, "Pct_of_Total"],
    pct_genes = bp_df[i, "Pct_of_Total_genes"]
  )
}

# Assign BP labels in same order as Python (sorted by height then node_id)
bp_labels_list <- vector("list", length(bp_indices))
for (j in seq_along(bp_indices)) {
  bp_id <- sprintf("BP%02d", j)
  info  <- bp_lookup[[ bp_id ]]
  bp_labels_list[[j]] <- list(
    merge_row  = bp_indices[j],
    height     = node_heights[ bp_indices[j] ],
    bp_id      = bp_id,
    pct_mods   = if (!is.null(info)) info$pct_mods  else NA,
    pct_genes  = if (!is.null(info)) info$pct_genes else NA
  )
}

# ---- 4. Compute x-positions for each internal node ----
# hclust leaf order: hc$order gives left-to-right positions (1-indexed)
leaf_pos <- setNames(seq_along(hc$order) - 0.5, hc$order)
# leaf_pos[leaf_index] = x position in [0, n]

node_x_pos <- numeric(n - 1)   # for internal nodes (merge rows)

# Iterative x-position computation (avoids stack overflow from recursion).
# hclust guarantees merge[i,] only references rows with index < i for internal
# nodes, so processing rows in order is safe.
for (i in seq_len(n - 1)) {
  left  <- merge[i, 1]
  right <- merge[i, 2]
  xl <- if (left  > 0) leaf_pos[left]  else node_x_pos[-left]
  xr <- if (right > 0) leaf_pos[right] else node_x_pos[-right]
  node_x_pos[i] <- (xl + xr) / 2
}

# ---- 5. Draw heatmap ----
# Reorder columns of eig_mat to match dendrogram leaf order
ordered_modules <- modules[ hc$order ]
eig_ordered     <- eig_mat[, ordered_modules]

# Color scale
vmax <- quantile(abs(eig_ordered), 0.98)
col_fun <- colorRamp2(
  c(-vmax, 0, vmax),
  c("blue", "white", "red")
)

ht <- Heatmap(
  eig_ordered,
  name            = "Eigengene",
  col             = col_fun,
  cluster_rows    = FALSE,
  cluster_columns = as.dendrogram(hc),   # use our hclust dendrogram
  column_dend_height = unit(6, "cm"),
  show_column_names  = TRUE,
  column_names_rot   = 90,
  column_names_gp    = gpar(fontsize = 7),
  row_labels         = rownames(eig_mat),
  row_names_gp       = gpar(fontsize = 9),
  row_title          = "Sample",
  column_title       = "Jorstad metamodules, consensusMean",
  column_title_gp    = gpar(fontsize = 15, fontface = "plain"),
  heatmap_legend_param = list(
    title            = "Eigengene value",
    direction        = "horizontal",
    title_position   = "topcenter",
    legend_width     = unit(4, "cm"),
    labels_gp        = gpar(fontsize = 9),
    title_gp         = gpar(fontsize = 11)
  )
)

# ---- 6. Draw and annotate branchpoint labels ----
pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/exploratory/jorstad_consensusMean_labeledWithClaude"), "dendrogram20_complexheatmap.pdf")#, 
#width = 4000, height = 2200
)

draw(ht,
     heatmap_legend_side = "bottom",
     padding = unit(c(8, 2, 2, 2), "mm"))   # bottom, left, top, right

# After draw(), add branchpoint labels onto the dendrogram viewport
# The column dendrogram is drawn in viewport "root/heatmap_Eigengene/column_dend"
# We use a decorate_* call to add grobs inside the dendrogram panel

decorate_column_dend("Eigengene", {

  # Current viewport dimensions
  vp_width  <- convertWidth(unit(1, "npc"),  "mm", valueOnly = TRUE)
  vp_height <- convertHeight(unit(1, "npc"), "mm", valueOnly = TRUE)

  # x in viewport: columns span 0..n in data space mapped to 0..1 npc
  # node_x_pos is in [0, n]; scale to npc
  x_scale <- 1 / n

  # y in viewport: dendrogram height spans 0..max_height mapped to 0..1 npc
  max_h   <- max(hc$height)
  y_scale <- 1 / max_h

  for (bp in bp_labels_list) {
    x_npc <- node_x_pos[ bp$merge_row ] * x_scale
    y_npc <- bp$height * y_scale

    label_text <- sprintf("%s\n%.1f%% (m)\n%.1f%% (g)",
                          bp$bp_id, bp$pct_mods, bp$pct_genes)

    grid.text(
      label      = label_text,
      x          = unit(x_npc, "npc"),
      y          = unit(y_npc, "npc"),
      just       = c("center", "bottom"),
      gp         = gpar(
        fontsize   = 7,
        fontface   = "bold",
        col        = rgb(0.55, 0, 0, 0.75)
      )
    )

    # Optional: background box
    grid.rect(
      x      = unit(x_npc, "npc"),
      y      = unit(y_npc, "npc"),
      width  = unit(1.8, "cm"),
      height = unit(0.9, "cm"),
      just   = c("center", "bottom"),
      gp     = gpar(fill = rgb(1, 1, 0.9, 0.5), col = NA)
    )

    # Re-draw text on top of box
    grid.text(
      label      = label_text,
      x          = unit(x_npc, "npc"),
      y          = unit(y_npc, "npc"),
      just       = c("center", "bottom"),
      gp         = gpar(
        fontsize   = 7,
        fontface   = "bold",
        col        = rgb(0.55, 0, 0, 0.75)
      )
    )
  }

  # Cut line
  grid.lines(
    x  = unit(c(0, 1), "npc"),
    y  = unit(rep(cut_height * y_scale, 2), "npc"),
    gp = gpar(col = "red", lty = 2, lwd = 2)
  )
})

dev.off()
message("Saved: dendrogram20_complexheatmap.pdf")

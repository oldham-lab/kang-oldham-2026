# ============================================================
# build_branchpoint_table()
#
# Computes a hierarchical clustering on the columns of a
# module eigengene matrix (Pearson correlation, complete
# linkage, distance = 1 - r) and returns a data frame of
# all branchpoints at or above a cut height, with their
# member modules and percentage of total.
#
# Optionally writes the table to a CSV file.
#
# Arguments:
#   mod_eig_file  : path to mod_eig.csv (rows = cell types,
#                   cols = modules, first col = row names)
#   cut_height    : minimum merge height to label as a
#                   branchpoint (default 0.3)
#   output_file   : if not NULL, path to write CSV output
#
# Returns:
#   A data frame with columns:
#     Branchpoint_ID      – BP01, BP02, …
#     Merge_Height_(1-r)  – height of the merge
#     N_Elements          – number of leaf modules
#     Pct_of_Total        – N_Elements / total modules * 100
#     Elements            – comma-separated module names
#
# Required packages: none (base R only)
#
# Example:
#   tbl <- build_branchpoint_table(
#     mod_eig_file = "mod_eig.csv",
#     cut_height   = 0.3,
#     output_file  = "branchpoint_table_modeig.csv"
#   )
# ============================================================

build_branchpoint_table <- function(mod_eig_file,
                                    cut_height  = 0.3,
                                    output_file = NULL) {

  # ---- 1. Load data ----
  eig     <- read.csv(mod_eig_file, header = TRUE, row.names = 1,
                      check.names = FALSE)
  modules <- colnames(eig)
  n       <- length(modules)

  # ---- 2. Pearson correlation -> distance -> complete linkage ----
  corr_mat <- cor(as.matrix(eig), method = "pearson")
  dist_mat <- as.dist(1 - corr_mat)
  hc       <- hclust(dist_mat, method = "complete")

  # ---- 3. Get all leaves descending from a given node (iterative) ----
  get_leaves <- function(target_node) {
    # hc$merge uses negative indices for original leaves,
    # positive for internal nodes (1-indexed)
    stack  <- target_node
    leaves <- integer(0)
    while (length(stack) > 0) {
      node  <- stack[length(stack)]
      stack <- stack[-length(stack)]
      left  <- hc$merge[node, 1]
      right <- hc$merge[node, 2]
      for (child in c(left, right)) {
        if (child < 0) {
          leaves <- c(leaves, -child)   # leaf: original index
        } else {
          stack <- c(stack, child)      # internal node: recurse
        }
      }
    }
    sort(leaves)
  }

  # ---- 4. Identify merge nodes at or above cut_height ----
  n_merges    <- nrow(hc$merge)
  bp_indices  <- which(hc$height >= cut_height)

  # Sort by ascending height (ties broken by node index) to match Python labeling
  bp_indices  <- bp_indices[order(hc$height[bp_indices], bp_indices)]

  # ---- 5. Build output table ----
  rows <- lapply(seq_along(bp_indices), function(j) {
    i         <- bp_indices[j]
    bp_id     <- sprintf("BP%02d", j)
    leaves    <- get_leaves(i)
    leaf_names <- modules[leaves]
    n_leaves  <- length(leaves)
    pct       <- round(n_leaves / n * 100, 2)
    data.frame(
      Branchpoint_ID         = bp_id,
      `Merge_Height_(1-r)`   = round(hc$height[i], 6),
      N_Elements             = n_leaves,
      Pct_of_Total           = pct,
      Elements               = paste(leaf_names, collapse = ", "),
      stringsAsFactors       = FALSE,
      check.names            = FALSE
    )
  })

  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL

  # ---- 6. Optionally write to CSV ----
  if (!is.null(output_file)) {
    write.csv(tbl, file = output_file, row.names = FALSE, quote = TRUE)
    message("Saved: ", output_file)
  }

  return(tbl)
}


# ---- Example usage (uncomment to run) ----
# tbl <- build_branchpoint_table(
#   mod_eig_file = "mod_eig.csv",
#   cut_height   = 0.3,
#   output_file  = "branchpoint_table_modeig.csv"
# )
# print(head(tbl))

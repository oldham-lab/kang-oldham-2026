# PC1-variance-specific definitions. Sources the shared grouping/helpers from
# ../mod_group_defs.R, then runs per-module correlation PCA on the bulk matrix and
# exposes both a per-gene and a per-module view of PC1 variance explained.
#
# Provides:
#   r2_tbl   gene -> (mod, r2): each gene's % variance explained by its module's
#            PC1 (= cor(gene, module-PC1)^2 * 100). Celltype-agnostic (bulk).
#   pve_tbl  per module: n_genes + pve (= mean per-gene r2 = PC1's overall % var
#            explained for that module).
#   build_grouped_r2(celltypes)  gene-level r2 replicated across cell types and
#            grouped Shared/Unshared/Neither by each cell type's membership.

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mod_group_defs.R"))

bulk_mat <- as.matrix(bulk_expr[, 3:ncol(bulk_expr)])   # genes x samples
rownames(bulk_mat) <- bulk_expr[[2]]                    # gene symbols (dups -> first)

# Per-gene variance explained by the module's first PC (correlation PCA over samples).
pca_one <- function(genes) {
  g <- intersect(genes, rownames(bulk_mat))
  if (length(g) < 2) return(NULL)
  M <- t(bulk_mat[g, , drop = FALSE])                   # samples x genes
  M <- M[, apply(M, 2, sd) > 0, drop = FALSE]           # drop constant genes
  if (ncol(M) < 2) return(NULL)
  pc1 <- prcomp(M, center = TRUE, scale. = TRUE)$x[, 1]
  tibble(gene = colnames(M), r2 = apply(M, 2, function(col) cor(col, pc1)^2) * 100)
}

mod_genes <- split(gene_tbl$gene, gene_tbl$mod)
r2_tbl <- bind_rows(lapply(these_mods, function(m) {
  t <- pca_one(mod_genes[[as.character(m)]]); if (!is.null(t)) t$mod <- m; t
}))
pve_tbl <- r2_tbl %>%
  group_by(mod) %>%
  summarise(n_genes = n(), pve = mean(r2), .groups = "drop")

# Gene-level r2 replicated across cell types and grouped by membership (bulk is
# celltype-agnostic, so the same gene r2 is re-grouped per cell type).
build_grouped_r2 <- function(celltypes) {
  r2_tbl %>%
    crossing(celltype = celltypes) %>%
    mutate(group = group_in_celltype(mod, celltype), series = "Bulk") %>%
    filter(!is.na(group))
}

# Mean-expression-specific definitions. Sources the shared grouping/helpers from
# ../mod_group_defs.R, then builds z-scored mean expression for the single-cell
# datasets (Liu, Gabitto) and bulk, plus build_grouped() over a set of cell types.
#
# Provides:
#   sc_z                single-cell z-scored mean expression (Liu, Gabitto)
#   bulk_z              bulk z-scored mean expression (one series, all samples)
#   build_grouped(celltypes)  grouped, z-scored long table

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/mod_group_expression/mod_group_defs.R"))

# ---------------------------------------------------------------------------
# Single-cell mean expression -> z-scored within dataset
# ---------------------------------------------------------------------------
sc_sources <- tribble(
  ~dataset,  ~dir,
  "Liu",     file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output"), REGION, "means"),
  "Gabitto", file.path(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output"), REGION, "means")
)

read_sc <- function(dataset, dir, cond) {
  df <- fread(file.path(dir, sprintf("genomewide_means_%s.csv", cond)), data.table = F)
  names(df)[1] <- "gene"
  df %>%
    pivot_longer(-gene, names_to = "celltype", values_to = "expr") %>%
    mutate(dataset = dataset, condition = cond)
}

sc_z <- sc_sources %>%
  crossing(condition = c("Con", "allAD")) %>%
  pmap_dfr(function(dataset, dir, condition) read_sc(dataset, dir, condition)) %>%
  mutate(expr_scale = log1p(expr)) %>%               # linear means -> log1p
  group_by(dataset) %>%
  mutate(z = as.numeric(scale(expr_scale))) %>%      # z within dataset
  ungroup()

# ---------------------------------------------------------------------------
# Bulk mean expression (single series, all samples; celltype-agnostic)
# ---------------------------------------------------------------------------
bulk_z <- tibble(gene = bulk_expr[[2]],
                 expr = rowMeans(bulk_expr[, -c(1, 2)])) %>%
  mutate(expr_scale = expr,                          # ComBat already log-normalized
         z = as.numeric(scale(expr_scale)))

# ---------------------------------------------------------------------------
# Assemble grouped, z-scored long table over a set of cell types.
# Bulk is celltype-agnostic, so it is replicated across the requested cell types
# and re-grouped by each cell type's membership.
# ---------------------------------------------------------------------------
build_grouped <- function(celltypes) {
  sc <- sc_z %>%
    filter(celltype %in% celltypes) %>%
    inner_join(gene_tbl, by = "gene") %>%                  # adds mod
    mutate(group = group_in_celltype(mod, celltype))

  bulk <- bulk_z %>%
    inner_join(gene_tbl, by = "gene") %>%
    crossing(celltype = celltypes) %>%
    mutate(dataset = "Bulk", condition = "Bulk",
           group = group_in_celltype(mod, celltype))

  bind_rows(sc, bulk) %>%
    filter(!is.na(group)) %>%
    mutate(series = ifelse(dataset == "Bulk", "Bulk", paste(dataset, condition)))
}

# Module-group PC1 VARIANCE analysis (fig_8)

Per-module **PC1 % variance explained** (bulk only), compared across the
Shared / Unshared / Neither module groups. Same cell-type-resolved grouping as
the sibling `../mean_expression/` analysis; grouping + helpers + raw `bulk_expr`
come from `../mod_group_defs.R`. See `../CONTEXT.md` for the container overview.

## Metric

For each module: take the bulk expression of its genes (genes × samples), run
**correlation PCA** over samples, and for each gene compute the % of its variance
explained by PC1 (`cor(gene, PC1)^2 * 100`, the per-gene **R²**). The per-module
**PVE** is the mean per-gene R² (= PC1's overall proportion of variance explained
with correlation PCA). Computed once in `pc1_variance_defs.R`, which exposes both
the per-gene (`r2_tbl`) and per-module (`pve_tbl`) views. Modules with < 2
measured/non-constant bulk genes get NA (none here; all 1016 resolved).

Two figures use these two views:
- **Main** (`pc1_variance_by_celltype`): unit = module, value = per-module PVE.
- **Featured** (`pc1_variance_featured`): unit = gene, value = per-gene R²; mirrors
  `../mean_expression/shared_module_featured_bulk` (each shared module's genes vs
  Unshared & Neither genes), one panel per shared triplet.

## Key decisions

- **Bulk only** (single-nucleus ignored, per request).
- **Cell-type-resolved + faceted** (mirrors the bulk mean-expression figure): each
  module's PVE is replicated across all 23 cell types and re-grouped by that cell
  type's membership. Consequence: a cell type's **Shared** group is tiny (1–7
  modules; the 18 shared triplets total), so a Shared bracket appears only where
  ≥ 3 shared modules exist (2 cell types). Unshared/Neither groups are large.
- Main figure unit = **module** (one point); featured figure unit = **gene**
  (per-gene R², so the featured module's genes form a distribution).
- Effect size = **Cliff's delta**; both figures carry pairwise Wilcoxon brackets
  (BH-adjusted p + stars), placed by `bracket_layout()` (shared defs) in headroom
  above each facet's whiskers and stacked, as in the mean-expression bulk figures.

## Inputs

- Grouping + `gene_tbl` + `expr_celltypes` + `bulk_expr`: from `../mod_group_defs.R`.
- Bulk matrix: `combined_FCX_final_1_1518_ComBat.csv` (genes × 1518 samples, ComBat-log).

## Files

| File | Description |
|---|---|
| `pc1_variance_defs.R` | Per-module correlation PCA: builds `r2_tbl` (per-gene R²), `pve_tbl` (per-module PVE), and `build_grouped_r2()`. Sources `../mod_group_defs.R`; sourced by both figure scripts. |
| `pc1_variance_figure.R` | Main figure: per-module PVE replicated across cell types, Wilcoxon + faceted bulk figure with brackets. |
| `pc1_variance_by_celltype.{png,pdf}` | Main figure (23 cell-type facets, Shared/Unshared/Neither). |
| `pc1_variance_featured_figure.R` | Featured figure: each shared triplet's module genes ("Module") vs Unshared & Neither genes, value = per-gene R². |
| `pc1_variance_featured.{png,pdf}` | 18-panel featured figure (one per shared triplet). |
| `pc1_variance_per_module.csv` | Per-module `mod`, `n_genes`, `pve`. |
| `pc1_variance_per_gene.csv` | Per-gene `gene`, `r2`, `mod`. |
| `pc1_variance_long.csv` | Main-figure tidy data (mod, n_genes, pve, celltype, group, series). |
| `pc1_variance_wilcoxon.csv` | Main-figure tests (`cliffs_delta`, `median_diff`, `p`, `p_adj_BH`). |
| `pc1_variance_featured_long.csv` | Featured tidy data (gene, r2, mod, celltype, group, plot_group, panel, …). |
| `pc1_variance_featured_wilcoxon.csv` | Featured per-triplet tests with effect sizes. |
| `pc1_variance_allct_figure.R` | "Over all mods", cell-type-agnostic: single pooled panel (no facet), per-module PVE grouped by `group_overall` (Shared/Unshared/Neither defined across any cell type; see `../mod_group_defs.R`). Each module contributes ONE value to ONE group (no cell-type replication). Wilcoxon brackets. |
| `pc1_variance_allct.{png,pdf}`, `pc1_variance_allct_long.csv`, `pc1_variance_allct_wilcoxon.csv` | Over-all-mods figure + backing data/tests. |
| `pc1_variance_featured_allct_figure.R` | "Per mod", cell-type-agnostic: mirrors `pc1_variance_featured_figure.R` but one panel per shared MODULE (12 unique modules, not 18 triplets); per-gene R² vs Unshared & Neither genes. Unshared/Neither are the same set in every panel. |
| `pc1_variance_featured_allct.{png,pdf}`, `pc1_variance_featured_allct_long.csv`, `pc1_variance_featured_allct_wilcoxon.csv` | Per-mod 12-panel figure + backing data/tests. |

Run order: each figure script sources `pc1_variance_defs.R` (which sources `../mod_group_defs.R`) automatically.

## Cell-type-agnostic ("over all cell types") figures

`pc1_variance_allct_figure.R` ("over all mods", unit = module) and
`pc1_variance_featured_allct_figure.R` ("per mod", unit = gene) collapse the
cell-type dimension that the main/featured figures facet on, using `group_overall`
(`../mod_group_defs.R`). Bulk only. Over-all-mods: Shared (n=12) > Unshared (n=84)
> Neither (n=920) modules, all three pairs significant. Per-mod: Module > Neither
per-gene R² in all 12 panels.

## Results

- All **1016 modules** resolved; **median PVE ≈ 70%**.
- **Main**: 21 Wilcoxon tests, 18 significant at BH < 0.05. Shared > Neither in 2/2
  strata where a Shared group is testable; Unshared > Neither broadly (DE modules
  are more cohesive / co-expressed in bulk than non-DE modules).
- **Featured**: 18 panels, 54 tests, 47 significant. Featured (shared) module's
  genes > Neither genes in per-gene PC1 R² in 18/18 panels.

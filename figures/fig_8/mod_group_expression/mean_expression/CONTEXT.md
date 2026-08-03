# Module-group MEAN EXPRESSION analysis (fig_8)

Compares **mean gene expression** across **module groups** defined by the overlap
of the AD and SCZ DE objects, in three expression datasets. Grouping + shared
helpers come from `../mod_group_defs.R`; see `../CONTEXT.md` for the container
overview and the sibling `../pc1_variance/` analysis.

## "Shared" definition (current)

Shared is defined per **(module, cell type, direction)** triplet: a triplet
present in BOTH the AD and SCZ DE tables. There are **18** such triplets (all
direction = down). Grouping is therefore **cell-type-resolved**:

Within a given cell type *c*:
- **Shared**   — module is a shared triplet in *c*
- **Unshared** — module is DE in *c* in only one disease (or opposite directions)
- **Neither**  — module (in the 1016 universe) is not DE in *c*

A module can be Shared in one cell type and Neither in another.

> An earlier version defined Shared at the module-ID level (32 modules shared by
> ID, ignoring cell type/direction). That was replaced by the triplet definition.

## Inputs

- AD / SCZ DE tables: `fig_7/dcopa_svg_map_path/{ad,scz}_dfc.csv` (cols: mod, Celltype, Direction, …)
- Module → gene map (`topmodposbc`): `…/SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv`
- 1016-module universe: size > 3 and bulk-cor support ≥ 2 (`bulk_cors_sigcount_bonf_1158.csv`), reproducing fig_3.
- Expression (per-gene means, region = **PFC**):
  - **Liu (MIT)**: `…/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/means/genomewide_means_{Con,allAD}.csv`
  - **Gabitto (SeAAD2024)**: `…/SEAAD2024_full_python_output/PFC/means/genomewide_means_{Con,allAD}.csv`
  - **Bulk**: `combined_FCX_final_1_1518_ComBat.csv` (fig_3 line 23); one series, mean over all 1518 samples.

## Key decisions

- Conditions: **Con + allAD** for the single-cell datasets.
- Cell type kept as a **facet** (not collapsed).
- **Z-score within each dataset**: single-cell means `log1p`'d first (linear, right-skewed);
  bulk is ComBat-log already; then z within dataset. Wilcoxon run on `expr_scale`
  (within-dataset z is a linear transform → same p).
- AD/SCZ cell-type names map to expression columns by case-insensitive match
  (only interneuron labels differ, e.g. `PVALB` → `Pvalb`).
- Bulk is celltype-agnostic: replicated across cell types and re-grouped by each
  cell type's membership.
- **Bulk and single-nucleus rendered as separate figures** (each with its own
  y-axis zoom), so bulk is no longer compressed on a shared z-axis.
- Main figure shows **all 23 cell types**; cell types with no shared triplet have
  only Unshared / Neither (no Shared box). The 5 with a shared triplet (L5 IT,
  L5/6 NP, Lamp5 Lhx6, Pvalb, Sncg) have all three groups.
- Effect size = **Cliff's delta** (= rank-biserial r, [-1, 1]); `median_diff` also reported.
- **Bulk figures carry pairwise Wilcoxon brackets** (`geom_signif`): stars + the
  BH-adjusted p for each group pair. Brackets are placed by `bracket_layout()`
  (shared defs) in headroom **above** each facet's whiskers and stacked, and the
  bulk z is winsorized at the 97th pct for display so boxes stay legible without
  whiskers colliding with the brackets.
- **Single-nucleus stats** are shown as a **companion effect-size heatmap**
  (`delta_heatmap`): rows = cell type (or triplet), x = series, facet = comparison;
  tile = Cliff's delta, text = BH-adj p stars. (Brackets don't work on the sn
  boxplots — 4 dodged series per group make them ambiguous.)

## Files

| File | Description |
|---|---|
| `../mod_group_defs.R` | Shared definitions (parent folder): 1016 universe, gene→module, 18 shared triplets, `group_in_celltype()`, raw `bulk_expr`, plotting helpers. |
| `mean_expression_defs.R` | Mean-expression transforms: `sc_z`, `bulk_z`, `build_grouped()`. Sources `../mod_group_defs.R`; sourced by both figure scripts. |
| `mod_group_expression_figure.R` | Main figure + Wilcoxon (Shared/Unshared/Neither per dataset × condition × cell type), all 23 cell types; renders separate sn and bulk figures. |
| `mod_group_expression_by_celltype_{sn,bulk}.{png,pdf}` | Main figure, single-nucleus and bulk versions (bulk carries Wilcoxon brackets). |
| `mod_group_expression_by_celltype_sn_stats.{png,pdf}` | Companion Cliff's-delta + stars heatmap for the sn main figure. |
| `mod_group_expression_long_{sn,bulk}.csv` | Tidy backing data (gene, mod, group, dataset, condition, celltype, expr, expr_scale, z, series), split sn vs bulk. |
| `mod_group_expression_wilcoxon_{sn,bulk}.csv` | Tests with `cliffs_delta`, `median_diff`, `p`, `p_adj_BH`, split sn vs bulk. |
| `shared_module_featured_figure.R` | Per-triplet figure: each shared (mod, ct, dir) featured as "Module" vs Unshared & Neither genes in that cell type; renders separate sn and bulk figures. |
| `shared_module_featured_{sn,bulk}.{png,pdf}` | 18-panel figure (one per shared triplet), single-nucleus and bulk versions (bulk carries Wilcoxon brackets). |
| `shared_module_featured_sn_stats.{png,pdf}` | Companion Cliff's-delta + stars heatmap (rows = triplet) for the sn per-triplet figure. |
| `shared_module_featured_long_{sn,bulk}.csv` | Tidy backing data for the per-triplet figure, split sn vs bulk. |
| `shared_module_featured_wilcoxon_{sn,bulk}.csv` | Per-triplet tests with effect sizes, split sn vs bulk. |
| `mod_group_expression_allct_figure.R` | "Over all mods", cell-type-agnostic: one pooled panel (no cell-type facet), Shared/Unshared/Neither by `group_overall` (a module is Shared if a shared triplet in ANY cell type, Unshared if DE in any, else Neither). sc expression pooled over all 23 cell types (every gene×cell type value is a point); bulk one value per gene. Renders separate sn + bulk figures + sn stats heatmap. |
| `mod_group_expression_allct_{bulk,sn}.{png,pdf}` | "Over all mods" figure, bulk (Wilcoxon brackets) and single-nucleus (4 dodged series). |
| `mod_group_expression_allct_sn_stats.{png,pdf}` | Companion Cliff's-delta heatmap for the sn over-all-mods figure. |
| `mod_group_expression_allct_long_{sn,bulk}.csv`, `mod_group_expression_allct_wilcoxon_{sn,bulk}.csv` | Backing data + tests for the over-all-mods figure. |
| `shared_module_featured_allct_figure.R` | "Per mod", cell-type-agnostic: mirrors `shared_module_featured_figure.R` but one panel per shared MODULE (12 unique modules, not 18 triplets), expression pooled over all cell types; Unshared/Neither are the same set in every panel (only the featured module changes). Renders separate sn + bulk figures + sn stats heatmap. |
| `shared_module_featured_allct_{bulk,sn}.{png,pdf}` | "Per mod" 12-panel figure, bulk (brackets) and single-nucleus. |
| `shared_module_featured_allct_sn_stats.{png,pdf}` | Companion Cliff's-delta heatmap (rows = featured module) for the sn per-mod figure. |
| `shared_module_featured_allct_long_bulk.csv`, `shared_module_featured_allct_summary.csv`, `shared_module_featured_allct_wilcoxon_{sn,bulk}.csv` | Backing data for the per-mod figure. The full sn long table (~19M rows) is not written — it is the over-all-mods sn long data re-grouped per featured module; a per-(panel, group, series) summary is saved instead. |

Run order: either figure script sources `mean_expression_defs.R` (which sources `../mod_group_defs.R`) automatically.

## Cell-type-agnostic ("over all cell types") figures

`mod_group_expression_allct_figure.R` ("over all mods") and
`shared_module_featured_allct_figure.R` ("per mod") collapse the cell-type
dimension that the per-cell-type figures facet on: groups use `group_overall`
(`../mod_group_defs.R`), and single-nucleus expression is **pooled** over all 23
cell types (every gene×cell type value is a point). Pooling inflates n, so the
Wilcoxon p-values are tiny — **Cliff's delta is the meaningful effect size**
(Shared/Module > Neither ≈ 0.5, all 12 per-mod panels Module > Neither in bulk).
The Cliff's-delta computation coerces `n1*n2` to numeric to avoid integer overflow
at the large pooled n.

## Results

- **Main**: 155 tests (all 23 cell types), 138 significant at BH < 0.05. Shared >
  Neither in 25/25 strata (Cliff's δ ≈ 0.6). Shared vs Unshared mostly n.s. — both
  are DE modules; the gap is vs Neither.
- **Per-triplet**: 18 panels, 270 tests, 216 significant. Featured module >
  Neither in 90/90 strata (every shared module, every dataset/condition).

## Open items

- `ad_scz_comparison.R` (this folder; writes `../ad_scz_mod_celltype_direction_intersect.csv`)
  Fisher test is still module-level (32 / 1016); not redone at the triplet level.

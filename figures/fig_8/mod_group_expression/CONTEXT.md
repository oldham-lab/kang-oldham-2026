# Module-group analyses (fig_8)

Container for analyses that compare a per-module / per-gene quantity across
**module groups** defined by the overlap of the AD and SCZ DE objects
(Shared / Unshared / Neither). The grouping, module universe, gene→module map,
plotting helpers, and the raw bulk expression matrix are shared; each analysis
lives in its own subfolder with an analysis-specific `*_defs.R` that sources the
shared `mod_group_defs.R`.

## Shared files (this folder)

| File | Description |
|---|---|
| `mod_group_defs.R` | Shared definitions: 1016-module universe, `gene_tbl`, 18 shared (mod, ct, dir) triplets, `group_in_celltype()`, `group_overall()` (cell-type-agnostic group: Shared if a shared triplet in any cell type, Unshared if DE in any, else Neither — 12/84/920 modules), `expr_celltypes`, raw `bulk_expr`, and plotting helpers (`sig_stars`, `signif_brackets`, `bracket_layout`, `delta_heatmap`). Sourced by each subfolder's `*_defs.R`. |
| `ad_scz_comparison.R` | Precursor: builds the AD/SCZ overlap, the 18-triplet intersect (`../ad_scz_mod_celltype_direction_intersect.csv`), and the module-level Fisher test (32 / 1016). |

## "Shared" definition

Per **(module, cell type, direction)** triplet present in BOTH the AD and SCZ DE
tables → **18 triplets** (all direction = down). Grouping is cell-type-resolved:
within a cell type, a module is **Shared** (a shared triplet there), **Unshared**
(DE in only one disease), or **Neither** (in the 1016 universe, not DE).

## Subfolders

| Subfolder | Analysis |
|---|---|
| `mean_expression/` | Mean gene expression per module group, in Liu + Gabitto (single-nucleus, Con & allAD) and bulk. Two figure families: main (Shared/Unshared/Neither, all 23 cell types) and per-triplet "featured" (18 panels). See `mean_expression/CONTEXT.md`. |
| `pc1_variance/` | Per-module PC1 % variance explained (bulk only): main (per-module PVE, faceted by cell type) + featured (per-gene R², 18 panels), mirroring the mean-expression figures. See `pc1_variance/CONTEXT.md`. |

## Shared display conventions (apply to both analyses)

- **Bulk vs single-nucleus split into separate figures**, each with its own
  y-axis, and **all 23 cell types** shown (cell types with no shared triplet get
  only Unshared / Neither).
- **Stats on bulk figures = in-plot Wilcoxon brackets** (`geom_signif`): stars +
  BH-adjusted p per group pair. `bracket_layout()` places them in headroom
  **above** each facet's boxplot whiskers and stacks them with enough room for the
  two-line labels (no overlap with boxes/whiskers or each other). Mean bulk z is
  winsorized at the 97th pct for display so boxes stay legible (PC1 R² is already
  bounded, so it is not winsorized).
- **Stats on single-nucleus figures = companion `delta_heatmap`** (Cliff's delta
  tiles + stars), because 4 dodged series per group make in-plot brackets
  ambiguous. (pc1_variance is bulk-only, so it has no sn figures/heatmaps.)
- Effect size = **Cliff's delta** throughout; Wilcoxon run on the analysis value
  (expr_scale / pve / per-gene r2).
- Figures use enlarged fonts (main `base_size` 17, featured 15, heatmap 16) on
  taller canvases (main 16×16, featured 13×21) so the dense bracket facets stay
  legible.

## Cell-type-agnostic ("over all cell types") figures

Both analyses gained two figures that collapse the cell-type facet (groups via
`group_overall`): an **"over all mods"** figure (all modules pooled into
Shared/Unshared/Neither in one panel) and a **"per mod"** figure (one facet per
shared MODULE — 12 unique modules, vs the 18 triplets of the per-cell-type
"featured" figures). For mean_expression the single-nucleus values are pooled over
all 23 cell types, so n is huge and Cliff's delta (not p) is the meaningful effect
size; pc1_variance is bulk only. Files: `*_allct*` in each subfolder.

## Session history (latest first)

- Added the cell-type-agnostic `*_allct` figures (over-all-mods + per-mod) to both
  analyses; added `group_overall()` to `mod_group_defs.R`.
- Increased font sizes across all figures; widened bracket spacing + enlarged
  canvases so the 3-group bracket facets stay clean.
- Added `pc1_variance/` featured figure (per-gene PC1 R², mirroring the mean
  featured bulk figure); factored its PCA into `pc1_variance_defs.R`.
- Fixed Wilcoxon bracket/whisker/text overlap on all bulk figures via
  `bracket_layout()` + winsorizing.
- Built the `pc1_variance/` analysis (per-module PC1 variance explained, bulk).
- Reorganised into this container: shared `mod_group_defs.R` + `mean_expression/`
  and `pc1_variance/` subfolders, each with its own `*_defs.R`.
- mean_expression: split bulk/sn into separate figures, extended to all 23 cell
  types, added bulk brackets + sn companion heatmap, split backing CSVs `_sn`/`_bulk`.

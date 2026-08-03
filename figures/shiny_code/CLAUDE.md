# CoPA Shiny App — Development Context

## Key Files

- **Main app**: `/home/gugene/ShinyApps/CoPA/app.R`
- **www root**: `/home/gugene/ShinyApps/CoPA/www/`

## App Structure (tabs)

1. **CoPA** — module projection plots (5 bslib cards in a CSS grid)
2. **dCoPA** — differential CoPA dot plots + SVG snapshots
3. **Core GBmap** — GBmap expression bar plots (plotly, two-column bslib cards)
4. **Gene Projection** — single-gene expression across datasets (ggplot, two-column bslib cards)
5. **About**
6. **How to Cite**

## Data / Resource Paths

All data files served through Shiny must live under `www/` (or be registered via `addResourcePath`):

- `www/topmodposbc_table.csv` — module KME table
- `www/bulk_cors_sigcount_bonf_1158.csv` — significance filter
- `www/unmerged_modules.qs` — seed modules
- `www/gsea_plot_dfs.qs` — GSEA data
- `www/CoPA_files/` — bulk_cor, mod_eig_var_explained, expr_line, projs, proj_cors
- `www/core_gbmap/` — annotation_level_{1-4}/ (means/var/n CSVs) + `calc_module_stats_core_gbmap.R`
- `www/dcopa/dot_data/` — 6 CSVs (ad_dfc, ad_dfc_admods, ad_mtg, ad_mtg_admods, scz_dfc, scz_dfc_sczmods)
- `www/dcopa/svg_map/` — same 6 CSV names
- **dCoPA SVG snapshots**: `www/dCoPA_snapshots/<context>/` — PNGs, served via the `svgs` alias. 6 app contexts: `AllADVsCon_DFC`, `AllADVsCon_DFC_ROSMAP`, `AllADVsCon_MTG`, `AllADVsCon_MTG_ROSMAP`, `SCZvsCon_DFC_bulkmegaset`, `SCZvsCon_DFC_bulkmegaset_Brainseq`. Pushed from the Fig 6 source SVGs (`Code_for_figures/fig_6/dCoPA_snapshots/`) by `Code_for_figures/shiny_code/svg_to_png.sh`, which copies only app-referenced modules (per `www/dcopa/svg_map/`) then rasterizes SVG→PNG at 144 DPI in place.
- `www/gene_projection/{dataset}/PFC/` — means/var/n **CSVs** for SEAAD2024, MIT, Morabito2021, LeinA9 (source only; not loaded by the app — see fst store below)

### On-disk fst store (what the app actually loads)

To keep RAM/startup low (and stop shinyapps.io OOM crashes), the big gene×cell-type
matrices and the bulk-correlation list are **not** held in RAM. They live on disk as
`fst` and are queried per-request (`.fst_rows()` / `.bulk_module()` in app.R). Only tiny
metadata is loaded at startup.

- `www/fst/gp/<ds>_{means,var,pct}.fst` — gene-projection matrices
- `www/fst/gb/annotation_level_<n>_{means,var}.fst` — core GBmap matrices
- `www/fst/bulk_cor.fst` — all modules' bulk correlations, stacked + row-indexed
- `www/fst/store_meta.qs` — gene lists, cell types, per-subclass n & global means,
  native-value quantiles, bulk row index, and relative fst paths

**Build step:** after the source CSVs/qs above are assembled, run
`Code_for_figures/shiny_code/build_fst_store.R`. It converts them to the fst store and
validates. Once built, the source `gene_projection/**/*.csv`, `core_gbmap/**/*.csv`, and
`CoPA_files/bulk_cor/bulk_cor_list.qs` are redundant for the deployed app (keep in
staging for rebuilds; drop from the bundle to keep it small).

Resource aliases registered:
```r
addResourcePath("copa_www", "/home/gugene/ShinyApps/CoPA/www")
addResourcePath("svgs",     file.path("/home/gugene/ShinyApps/CoPA/www", "dCoPA_snapshots"))
```

## Key Technical Decisions

### bitmapType = "cairo"
`options(bitmapType = "cairo")` **must be set** at startup. The server is headless (no X11 display). Without Cairo, plotly's internal `png()` calls fail with `unable to open connection to X11 display ''`. This is unrelated to font rendering.

### showtext — removed
`library(showtext)` and `showtext_auto()` were removed. showtext is not needed: fonts in ggplot `element_text(family = "sans")` render fine through Cairo without it. showtext caused a DPI mismatch (showtext default 96 DPI vs Shiny's 72 DPI) that made fonts appear incorrectly sized.

### bindCache and font sizes
`bindCache()` stores the **ggplot object** (not the rendered PNG) in an in-memory cache. If the app is restarted between removing showtext/cairo and comparing fonts, the cache is cleared and fonts render correctly. Stale cache entries from old rendering settings can persist across code changes until R is restarted.

**If fonts look wrong after changing rendering settings: restart R to clear the in-memory cache.**

### renderPlot width/height
All CoPA tab `renderPlot` calls explicitly read both width and height from `session$clientData`:
```r
width  = function() { w <- session$clientData$output_out_expr_width;  if (is.null(w) || w <= 0) 400L else w },
height = function() { h <- session$clientData$output_out_expr_height; if (is.null(h) || h <= 0) 400L else h }
```
This ensures consistent dimension handling. Gene Projection plots (tabs 4) specify only height (not width) via `renderPlot(..., height = function() {...})`.

### Layout — CoPA tab (tab 1)
CSS grid inside `mainPanel(width = 10)`:
```
height: calc(100vh - 234px); min-height: 600px; min-width: 1067px;
grid-template-columns: 1fr 1fr;
grid-template-rows: repeat(12, 1fr);
```
Cards use bslib `card()` + `card_header()` + `card_body(style = "flex:1; min-height:0; padding:0; overflow:hidden;")` with `plotOutput(height="100%", width="100%")`.

The `234px` offset accounts for: sticky nav strip (~66px) + copa_content margin-top (36px) + title + gene oneliner + grid padding + fixed footer (80px).

### Layout — Tabs 3 & 4 (Core GBmap, Gene Projection)
Two-column layout using `fluidRow(column(6, ...), column(6, ...))`. Each column contains a bslib card:
```r
card(
  card_header("log UMI counts", style = "text-align:center; font-weight:normal; font-size:14px; padding-top:3px; padding-bottom:2px; background:#e2e8f0; color:#4a5568;"),
  card_body(style = "flex:1; min-height:0; padding:8px; overflow:hidden;",
    withSpinner(plotOutput(..., height = "calc(100vh - 245px)"))
  )
)
```

### Download button — Gene Projection (tab 4)
Button starts disabled:
```r
actionButton("gp_screenshot", "Download as PDF", icon = icon("file-pdf"),
             style = "width: 100%; font-size: 12px;", disabled = "disabled")
```
Server observer toggles state:
```r
observe({
  if (length(input$gp_genes) > 0) enable("gp_screenshot") else disable("gp_screenshot")
})
```

### Download button — Core GBmap (tab 3)
Same pattern — starts disabled, enabled by server observer when `input$kme_file` and `input$module` are non-null (existing code at server lines ~2075–2083).

## CSS Notes

- Global reset: `*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }` — Bootstrap class selectors (specificity > 0) still win over this.
- Bootstrap column gutters removed for CoPA mainPanel: `.tab-pane > div > .row > .col-sm-10 { padding-left: 0; padding-right: 0; }`
- Nav strip: sticky, logo injected via JS (`copa_www/copa_header_image.png`), 56px height
- Fixed footer: 80px, `body { padding-bottom: 80px }` compensates
- `body { min-width: 1560px; overflow-x: auto; }`

# Tables S6 & S7 — v2 (current)

Pseudobulk modeling of gene expression as a function of cell-type abundance.
- **Table S6** — Jorstad et al. (`s6_jorstad.xlsx`)
- **Table S7** — SEAAD 2024 / Gabitto et al. (`s7_seaad2024.xlsx`)

Each workbook has a `Legend` tab plus one tab per Donor × Platform × Region ×
Celltype, holding per-gene `mean_expr_pcntile`, `adj_r2`, `rmse`.

## What changed from v1
Table S7's SEAAD data is from **Brodmann area A9 = DLPFC (dorsolateral prefrontal
cortex)**, but the v1 workbook mislabeled the region as **"MTG"** (a hard-coded
`region = "MTG"` in the generator; the loaded metadata is `SEAAD_A9_RNAseq…`).
v2 corrects the label to **DLPFC**. Only the label changed — all `adj_r2`/`rmse`
values are identical to v1.

## Files
| File | Description |
|---|---|
| `s6.csv`, `s7.csv` | Backing per-gene modeling results (gitignored; large). S7 region = DLPFC. |
| `s6_jorstad.xlsx`, `s7_seaad2024.xlsx` | The formatted tables (rendered from the CSVs). |
| `build_table_s6_7.R` | Rebuilds both `.xlsx` from the CSVs — **no raw data / model rerun needed**. |
| `table_6_7_v2.1.R` | Full generator: reruns the pseudobulk modeling from raw snRNA-seq to regenerate the CSVs + xlsx. `save_dir` points here. |

## Rebuild
```r
Rscript build_table_s6_7.R      # CSV -> xlsx (fast)
# or, to regenerate the CSVs from raw data first:
Rscript table_6_7_v2.1.R
```

## v1 (`../v1/`)
Archived pre-fix outputs: `s6_jorstad.xlsx`, `s7_seaad2024.xlsx` (S7 still labeled
"MTG"), and the older generator `table_6_7_v2.0.R`.

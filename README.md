# kang-oldham-2026

Analysis code for **Kang, Oldham, et al.** (working title). Reproduces the paper's
figures and tables from a batch-corrected bulk human-cortex RNA-seq megaset and
its projection onto single-nucleus datasets.

> Placeholder repo name. Method code lives in the companion package
> [**CoPA**](https://github.com/oldham-lab/COPA/tree/dev) (Co-expression Projection
> Analysis); this repo *uses* CoPA rather than re-implementing it.

## What it does

1. **Preprocessing** — assemble the ComBat-corrected bulk FCX megaset from ~7 cohorts
   (`preprocessing/`).
2. **Analyses** — derive consensus co-expression modules with CoPA and project them
   onto single-nucleus data (AD, SCZ vs control); real-vs-random / dCoPA significance
   (`analyses/`).
3. **Figures & tables** — the paper's figures and supplementary tables
   (`figures/`).

See `docs/copa-pipeline-flowchart.html` and `docs/consensus-analysis-WORKFLOW.md`
for the full run-order and data-flow.

## Layout

| Path | Contents |
|------|----------|
| `preprocessing/` | Stage 1: sample collection → quantification → QC/ComBat → megaset |
| `analyses/` | Stage 3–4: per-dataset CoPA application + random/dCoPA analysis |
| `figures/` | `fig_1`…`fig_8`, `fig_s*`, `table_s*` — one directory per figure/table |
| `docs/` | Pipeline flowchart, workflow write-up |

### Figure/table versions

Most figure and table directories are versioned (`v1`, `v2`, …). This repo keeps
only the versions the published figures actually need: **the latest version of
each figure, plus any earlier version that a kept file still reads from.** That
second clause is load-bearing — several latest versions read regenerated panels
or base workbooks out of an earlier one (`fig_3/v4`→`v3`, `fig_s11/v3`→`v2`,
`table_s1/v7`→`v6`, `table_s13/v3`→`v2`), and `fig_7/v3`, `v4`, `v7.1` and `v7.2`
are read by *other* figures and tables (`fig_8`, `fig_s11`, `table_s11`,
`table_s14`, `table_s15`). Superseded versions and exploratory scripts were
removed, as were the exploratory and supplementary directories
(`fig_2_related_analyses`, `fig_7_sup`, `fig_AD`, `fig_SCZ_AD`) and the Shiny
app pipeline (`shiny_code`), which lives with the app; the full history remains in the upstream
`oldham-lab/Consensus-analysis` repo.

## Requirements

- R (>= 4.1) with the **CoPA** package installed
  (`remotes::install_github("oldham-lab/COPA@dev")`). The `dev` ref is required:
  the repo's default branch holds the standalone engine scripts, not the package.
- Python (>= 3.9) for the single-nucleus processing bundled in CoPA (`inst/python`).
- See each figure directory for figure-specific dependencies.

## Configuration

Paths are centralized. Copy `.Renviron.example` to `.Renviron` and set `DATA_DIR`
and `REPO_DIR`, then `source("config.R")` at the top of a session. (Path
de-hardcoding across scripts is in progress — see the reorg notes.)

## Data availability

Datasets are **not** distributed with this repo (controlled-access AMP-AD/MSBB/
PsychENCODE cohorts via Synapse; public BrainSpan/SEA-AD). Cohorts, accessions,
and access models are documented in [`DATA.md`](DATA.md).

## Citation

If you use this code, please cite the paper and this software release — see
[`CITATION.cff`](CITATION.cff). The method package is
[CoPA](https://github.com/oldham-lab/COPA/tree/dev).

## License

GPL-3 (matches the CoPA dependency). See `LICENSE`.

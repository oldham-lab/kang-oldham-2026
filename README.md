# kang-oldham-2026

Analysis code for **Kang, Oldham, et al.** (working title). Reproduces the paper's
figures and tables from a batch-corrected bulk human-cortex RNA-seq megaset and
its projection onto single-nucleus datasets.

> Placeholder repo name. Method code lives in the companion package
> [**CoPA**](https://github.com/oldham-lab/COPA) (Co-expression Projection
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
| `_reorg_archive/` | Superseded/historical code (to be split into a separate archive repo) |

## Requirements

- R (>= 4.1) with the **CoPA** package installed (`remotes::install_github("oldham-lab/COPA")`).
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
[`CITATION.cff`](CITATION.cff). The method package is [CoPA](https://github.com/oldham-lab/CoPA).

## License

GPL-3 (matches the CoPA dependency). See `LICENSE`.

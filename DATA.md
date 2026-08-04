# Data availability

This repository ships **code, not data**. The datasets below are obtained from
their original providers and placed under the roots configured in `config.R` /
`.Renviron` (`DATA_DIR`, `MEGASET_DIR`, `DATA_OTHER_DIR`, `SHARED_DIR`, …). Some
are **controlled-access** and require an application/DUA with the provider.

> ⚠️ Accession IDs below are drawn from the download scripts in
> `preprocessing/1-sample_collection/`. **Confirm and complete this list before
> publication** (exact Synapse/SRA/dbGaP/AnVIL IDs, versions, and access model).

## Bulk RNA-seq cohorts (megaset — Stage 1 preprocessing)

| Cohort | Source | Access | Notes |
|--------|--------|--------|-------|
| BrainSeq | Synapse | controlled (AMP-AD/PsychENCODE) | + SCZ subset (`process_Brainseq_SCZ.R`) |
| GTEx | AnVIL (gen3) | controlled (dbGaP `phs000424`?) | `download_GTEx_from_anvil.sh` — confirm phs |
| NABEC | SRA / dbGaP | controlled | `phs001300`; `download_NABEC_from_SRA.sh` |
| ROSMAP | Synapse (AMP-AD) | controlled | `download_AMPAD_from_synapse` |
| MSBB | Synapse (AMP-AD) | controlled | frontal pole + IFG (`sif_MSBB_fp`, `sif_MSBB_ifg`) |
| CMC / CMC-HBCC | Synapse (PsychENCODE) | controlled | `download_CMC_from_synapse.py` |
| BrainGVEX | Synapse (PsychENCODE) | controlled | `download_psychEncode_from_synapse.py` |
| BrainSpan | SRA | public | `phs000755`?; `download_Brainspan_from_SRA.sh` |
| Li/Sestan 2018 | Synapse | controlled | `download_Li_Sestan_2018_from_synapse.py` |
| UCLA_ASD, Yale_ASD | (sample info only) | — | `sif_UCLA_ASD.csv`, `sif_Yale_ASD.csv` |

Sample-info tables (control-only, Age ≥ 18; SCZ subset) are tracked in
`preprocessing/sample_info_files/`. Full Synapse ID set referenced by the
download scripts: `syn3157322`, `syn3270014`, `syn3275211`, `syn3346807`,
`syn3354385`, `syn3817650`, `syn4566330`, `syn10476936`, `syn11038293`,
`syn11227080`, `syn11384571/589`, `syn11958660`, `syn12104376/381/384`,
`syn12299752`, `syn15672826`, `syn18134196`, `syn18358379`, `syn22331712`,
`syn22416880`, `syn24173489`, `syn35641127` (verify/curate).

## Single-nucleus datasets (Stage 3 projection)

| Dataset | Disease | Source | Access |
|---------|---------|--------|--------|
| SEA-AD 2024 (DFC/MTG) | AD vs control | Allen / AWS S3 | public |
| MIT Multiome (Gabitto labels) | AD vs control | Synapse (AMP-AD) | controlled |
| brainSCOPE (CMC, SZBDMulti-Seq) | SCZ vs control | brainscope.gersteinlab.org / PsychENCODE | see provider |

## Large reference files (not shipped)

- **MSigDB XML** for Broad GSEA — download from <https://www.gsea-msigdb.org/>;
  point `MSIGDB_XML` at it (used by CoPA and by `figures/fig_7/gsea_func*.R`).
- **Cell-type AD-vs-control DE results** (`fig_7`, ~140 MB of `.qs`/`.RDS`:
  `*_ADvsCon_by_celltype.{qs,RDS}` for MIT/SEA-AD × DFC/MTG). **Regenerate on
  demand — not tracked, not distributed.** Produced by
  `figures/fig_7/v4/full_DE_pipeline_ADvsCon.R` (writes per-dataset DE objects),
  then collated by `figures/fig_7/v4/DE/concat_DE_results.R`; downstream fig_7
  panels read them from `DE_DIR`. Set `DE_DIR` to the output location (default:
  `<REPO_DIR>/figures/fig_7/DE_old`, which is gitignored). Requires the
  single-nucleus inputs (SEA-AD, MIT Multiome) under `DATA_DIR`.
  (NB: `figures/fig_1/full_de_pipeline_from_scratch.R` is a *different* analysis —
  cell-type-vs-all marker DE — not the fig_7 AD-vs-control objects.)

## External analysis tools sourced by path

Three lab tools are `source()`d by absolute path rather than installed as
packages. Each has an env var (see `.Renviron.example`); clone or copy them and
point the var at your checkout:

| Var | Tool | Used by |
|---|---|---|
| `FINDMODULES_DIR` | FindModules (`FindModules/R/FindModules.R`) | `fig_4`, `fig_5`, `fig_s8`, fig_2 related analyses |
| `GSEA_GENERIC_DIR` | `GSEAfxsV3.r`, `GSEAfxsV3_nonpar_temp.r` | fig_2 related analyses (metacell/imputation FindModules runs) |
| `SAMPLENETWORK_DIR` | `SampleNetwork_1.08.r` | `preprocessing/3-process_counts/` |

> **⚠️ `SAMPLENETWORK_DIR` needs author confirmation.** The path these scripts
> originally referenced (`/home/gugene/code/SampleNetwork/`) no longer exists, and
> **two non-identical copies** of `SampleNetwork_1.08.r` survive on the analysis
> host: `/home/gugene/code/labcode_old/SampleNetwork/` and
> `/home/shared/code/SampleNetworks/`. The default points at the `labcode_old`
> copy, because that is the one `preprocessing/3-process_counts/1-cat_datasets.R`
> still referenced directly. Confirm which copy produced the published
> preprocessing before release, and pin it (ideally vendor it into this repo).

Smaller external inputs, each with its own var: `SCINRB_DIR` (scINRB imputation,
fig_2 related analyses only), `COMPAREMARKERS_DIR` (bulk fidelity table for the
fig_1 followup sensitivity analysis), `BBMAP_DIR` (`bbduk.sh` + adapter FASTA for
read trimming), `HOME_DATA_DIR` (two one-off lein MTG / ABI cell-count CSVs),
`PYTHON_BIN` (interpreter R invokes for SVG/PPTX rendering).

`SHINYAPP_DIR` points at the CoPA Shiny app checkout and is only used by
`figures/shiny_code/` snapshot generation — not by any figure.

## Tracked reproduction inputs (small)

- `analyses/bulk_module_significance/bulk_cors_sigcount_bonf_*.csv` — module
  significance pre-filter used across figures (relocated from the archived
  `fig_3_old/bulk_cor_significance_analysis/`; see that folder's README.md).
- `figures/fig_8/v1/scz_*.csv` + `scz_genes.txt` — curated SCZ gene database for
  fig_8, product of a non-deterministic LLM-agent literature-review pipeline
  (not cheaply regenerable). Tracked via a `.gitignore` exception; see
  `figures/fig_8/v1/README_scz_db.md`.

## Intermediate products

Module derivation, kME tables, projection indices, and figure outputs are
regenerated by the pipeline (see `docs/copa-pipeline-flowchart.html`) and are
**not** tracked (`.gitignore`). A subset of small analysis inputs (e.g. SCZ gene
databases under `figures/fig_8/`) are inputs to reproduction — see repo TODOs
for their final location.

## Configuration

Copy `.Renviron.example` → `.Renviron`, set the data roots, then `source("config.R")`.

# Workflow Outline — Consensus-analysis

> Describes the overall analysis workflow, including code sourced from outside this repo.

## What the whole thing does (one paragraph)

This repo builds a **consensus co-expression module reference** from a large bulk-brain RNA-seq "megaset" (~1,500 samples across 7 cohorts), then uses a **projection method (CoPA)** to map those modules onto single-nucleus RNA-seq data — establishing what cell type each module corresponds to — and finally runs a **differential projection (dCoPA)** to find which modules/cell types are dysregulated in disease (Alzheimer's and schizophrenia). Single-nucleus pseudobulk + differential-expression reanalysis provides the cell-type ground truth that validates the modules. Everything downstream (figures, tables, a Shiny app) visualizes one of these four products.

> **Naming note:** "CoPA" = **Co**-expression **P**rojection **A**nalysis; "dCoPA" = the **d**ifferential version. The README's phrase *"greedy march projection analyses"* = greedy-march modules (the reference) fed into CoPA (the projection). A few labels are unverified: what **"ESSA"** stands for (it's the DE-functions repo/dir name), and exact final module counts (scripts show ~1,300–1,400 *unmerged*, figures cite ~1,150 *merged* — approximate).

---

## 0. External code dependencies (sourced, NOT in this repo)

The repo is not self-contained — it sources six external toolkits. All exist on disk:

| Tool | Path | Role |
|------|------|------|
| **FindModules** | `/home/shared/code/FindModules/FindModules094.R` | Oldham-lab WGCNA-style co-expression module detection (clustering of a bicor/TOM matrix → modules, eigengenes, kME) |
| **SampleNetwork** | `/home/shared/code/SampleNetworks/SampleNetwork_1.08.r` | Sample-level QC: bicor outlier detection + ComBat batch correction |
| **COPA** | `/home/gugene/code/git/CoPA/` (entry point `wrapper.R`, sourced 68×) | The projection engine — builds modules, kME tables, projection indices, null distributions, GSEA, plotting objects |
| **GSEA_generic** | `/home/gugene/code/git/GSEA_generic/GSEAfxsV3*.r` | Fisher's-exact gene-set enrichment |
| **Pseudobulk-from-SC-SN-data** | `/home/gugene/code/git/Pseudobulk-from-SC-SN-data/makeSyntheticDatasets_0.51.r` | Generates synthetic bulk samples by sampling single cells |
| **ggplot_theme_settings** | `/home/gugene/code/ggplot_theme_settings.R` | Shared plot theming |

> ⚠️ **Path caveat:** Many `source()` calls point at `…/Consensus-analysis/RNAseq_analyses/1-FindModules/…`, `2-Projection/…`, `6-SN-pseudobulk/…`, `9-ESSA_DE/…`. **That `RNAseq_analyses/` tree no longer exists** — it was reorganized into `Analyses/` (e.g. `Analyses/FindModules_greedy_march/`, `Analyses/Projection_analyses/`, `Analyses/SN_pseudobulking/`, `Analyses/DE_reanalysis/`). So a chunk of the sourced paths are stale references to the current repo's own functions under their old layout.

---

## Stage A — Preprocessing: build the bulk "megaset"
`Preprocessing/` (stages 1 → 4)

1. **Collect** (`1-sample_collection/`, `1.1-disease_dataset_download/`) — download bulk cohorts: **BrainSeq, GTEx, NABEC, ROSMAP, MSBB, CMC, BrainGVEX** (+ ASD sets), from **Synapse / SRA / AnVIL / AWS S3**; SEA-AD single-nucleus from S3.
2. **Quantify** (`2-quantify_counts/`) — per-cohort `.sh`: FastQC/MultiQC → BBDuk trim → **Salmon** (BrainSeq) or **Kallisto** (most others) against a GRCh38/GENCODE index.
3. **Process** (`3-process_counts/`) — per-dataset `process_*.R`: `tximport` transcript→gene, **protein-coding filtering** (`protein_coding_gene_IDs.csv`), then **SampleNetwork** QC + ComBat. `1-cat_datasets.R` concatenates the 7 cohorts on shared genes → **`expr_combined.csv`** (~1,518 samples); `2-calculate_indiv_cor_mat.R` computes per-dataset **bicor similarity matrices**.
4. **QC analyses** (`4-analyses/`) — genome-wide bicor distributions, PC/covariate checks, normalization-strategy comparisons.

**Output → downstream:** combined ComBat-corrected expression matrix + per-dataset bicor similarity matrices.

---

## Stage B — Consensus module discovery ("greedy march")
`Analyses/FindModules_greedy_march/` + external `FindModules094.R`

- **Greedy march** (`1-run_greedy_march_FM.R`, `fxns/.../greedy_march_code_updated_20240802.R`): iterative module detection. Cluster the megaset correlation matrix, cut at a high threshold to extract modules, **remove assigned genes, lower the threshold, repeat** until a random-derived "floor." Produces ~1,300+ **unmerged** modules. (`2-FM_individual.R` runs standard FindModules per dataset for comparison.)
- **Merge** (`5-FM_alternate_merge_strats_application.R`): merge modules by **variance-explained** of the eigengene → final **consensus modules**.
- **FindCliques** (`4-FindCliques.R`) and `3`/`6`: alternate seed-finding (maximal cliques) and exploratory merge/enrichment strategies.

**Output → downstream:** consensus module definitions (gene→module), eigengenes, kME tables.

---

## Stage C — Single-nucleus pseudobulk + DE (cell-type ground truth)
`Analyses/SN_pseudobulking/` + `Analyses/DE_reanalysis/` + external `makeSyntheticDatasets`

- **Pseudobulk** SN/SC datasets (**Lein/Jorstad 2023** multi-region, **Bakken** 10x-vs-SSv4, **Miller**, **Martinowich/Tran 2020**, **SEA-AD 2024**) by aggregating nuclei per **cell type × donor**. `makeSyntheticDatasets` additionally builds synthetic bulk-from-SN samples to test cell-mixing effects.
- **DE reanalysis** (`DE_reanalysis/`, "ESSA"): **edgeR glmLRT** pseudobulk DE, one-vs-rest per cell type, blocked on donor, across Lein regions and SEA-AD (control/AD strata) → **cell-type marker gene lists**.
- **Cross-checks** (`5`–`8`): do bulk modules contain SN cell-type markers? Does FindModules-on-pseudobulk recapitulate bulk modules?

**Output → downstream:** pseudobulk cell-type × donor expression + DE marker sets — the reference CoPA projects onto.

---

## Stage D — Projection (CoPA): map modules → cell types
external `COPA/wrapper.R` driving `Analyses/Projection_analyses/`

Pipeline inside `wrapper.R`:
1. (re)derive greedy-march modules on the bulk reference;
2. **kME / topmodposbc** tables (assign each gene to its best module);
3. **GSEA** annotation of modules (`gsea.R` + GSEA_generic);
4. **summary tables** of the SN target (mean/var per cell type × condition);
5. **projection indices** — module "activity" per cell type, in 3 normalizations including the **REI (Relative Expression Index)**;
6. **null model**: project 10,000 random gene-sets → **euclidean distance + empirical p-values**.

**Output → downstream:** per-module-per-cell-type projection scores + significance = the module↔cell-type map.

---

## Stage E — Differential projection (dCoPA): disease dysregulation
`COPA/dCoPA_compare.R`, `Analyses/Projection_analyses/Brainseq_SCZsamples_COPA_runcode.R`, `calculate_rand_euclidean_distances.R`

Run CoPA on **disease vs. control** strata (AD: SEA-AD, ROSMAP; SCZ: BrainSeq), take the **euclidean distance between conditions per module/cell type**, test against the random-projection null → which modules are dysregulated, in which cell types, in which direction. Cross-disease overlap (AD ∩ SCZ).

---

## Stage F — Figures, tables, Shiny
`Code_for_figures/` (ignore `*_old`, `*_deprecated`)

The figure sequence traces the narrative arc:

- **Validate** modules = cell types: `fig_1`, `fig_2`, `fig_s2–s5`, `table_s6_7` — pseudobulk DE & per-donor modeling.
- **Annotate** modules: `fig_3`, `fig_4`, `fig_5` (REI schematic + dendrograms), `fig_6` (all significant modules), `fig_s6–s8`, `table_s9/s10` (module membership).
- **Disease**: `fig_7` (AD dCoPA), `fig_8` (SCZ dCoPA), `fig_s9` + `table_s11` (permutation null), `table_s12` (dysregulated gene lists). The shared-signature exploration (`fig_AD`, `fig_SCZ_AD`) is upstream only.
- **Cohort summary**: `fig_s1`, `table_s1`, `table_s8`.
- **Delivery**: the interactive app (CoPA / dCoPA / gene-projection tabs, served from precomputed GSEA caches) ships with the app itself, not this repo.
- **Enrichment backbone** (cross-cutting): `gsea_func_optimized.R`, `generic_enrichment_function.r`.

---

## End-to-end data flow

```
Public cohorts (Synapse/SRA/AnVIL/S3)
   └─[A] trim → Salmon/Kallisto → tximport → protein-coding filter
         → SampleNetwork+ComBat → cat → MEGASET (~1518 samples) + bicor matrices
              └─[B] greedy march + FindModules → merge by VE → CONSENSUS MODULES
                    │
   SN/SC datasets ──┤[C] pseudobulk (cell type × donor) + edgeR DE → cell-type markers
   (Lein, SEA-AD,   │
    Bakken, …)      │
                    └─[D] CoPA (COPA/wrapper.R): project modules → SN
                          → kME, REI indices, GSEA, random-null p-values
                          → MODULE ↔ CELL-TYPE MAP
                                └─[E] dCoPA: disease vs control (AD: SEA-AD/ROSMAP; SCZ: BrainSeq)
                                      → dysregulated modules/cell types + p-values
                                            └─[F] figures, supp tables, Shiny app
```

---

## Confidence

Stages, ordering, tools, and data flow are well-supported by the code. Lower-confidence items: the "ESSA" acronym, exact module counts, and the stale `RNAseq_analyses/` source paths (the functions still exist, just relocated under `Analyses/`).

# Projection Analyses — Code Inventory

Inventory of all code files under `Analyses/Projection_analyses/`, their associated
outputs, and a grouping into new / old / related / other.

## Architecture note

None of these scripts contain the analysis logic — they are thin **"runcode" drivers**
that `source()` shared functions from `/home/gugene/code/git/CoPA/`:

- `wrapper.R`
- `project_random_and_calculate_euclidean_standalone.R`
- `calculate_euclidean_distances_v2_standalone.R`

**All data outputs are written to external paths** (`/mnt/bdata/...greedy_march_pipeline_output/...`
and `~/figure_sup_analyses/...`), *not* inside this repo folder. The only non-code files
physically present in the tree are 2 log files. This folder is therefore code-only and is
**not self-contained / reproducible** without those external paths and the `COPA/` library.

---

## 1. New code (not in `older_projection_code/` or `related_analyses/`)

Current, consolidated drivers. The three `AllComparisons` scripts share an identical
structure (loop over region × comparison × module-type, calling
`project_rand_and_calculate_euclidean`).

| Script | What it runs | Associated output (external) |
|---|---|---|
| `Brainseq_SCZsamples_COPA_runcode.R` | `COPA()` + `COPA_compare()` on BrainSeq SCZ bulk megaset | `…/greedy_march_pipeline_output/Brainseq_SCZ/` (module dir: `kme_tables/`, `modules/`, etc.) |
| `brainSCOPE_python/brainSCOPE_AllComparisons_runcode.R` | SCZ vs Con, regions CMC & SZBDMulti-Seq, mod_types `bulk_megaset` + `brainseq_scz` | `…/brainSCOPE/brainscope_means_SE_output/{region}/euclidean_distances/{label}_{mod_type}_output_table.csv` |
| `MIT_Multiome_gabittoLabels/MITMultiome_AllComparisons_GabLabels_runcode.R` | PFC & MTC × {allAD/earlyAD/lateAD/APOE} × {bulk_megaset, rosmap} = 16 runs | `…/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/{region}/euclidean_distances/` |
| `SEAAD2024_gabittoLabels/SEAAD2024_AllComparisons_GabLabels_runcode.R` | PFC & MTC × 4 comparisons × 2 mod_types = 16 runs | `…/SEAAD2024_full_python_output/{region}/euclidean_distances/` |

### New but deprecated (`MIT_Multiome_gabittoLabels/deprecated/`)

Superseded by the single `MITMultiome_AllComparisons_GabLabels_runcode.R` above.

| Script | Output |
|---|---|
| `MITMultiome_MTC_AllADVsCon_GabLabels_runcode.R` | `…/MIT_ADMultiome_AllADVsCon_MTC_GabittoLabels_test/` |
| `MITMultiome_MTC_ROSMAP_AllADVsCon_GabLabels_runcode.R` | `…/MIT_ADMultiome_AllADVsCon_MTC_ROSMAP_GabittoLabels/` |
| `MITMultiome_PFC_AllADVsCon_GabLabels_runcode.R` | `…/MIT_ADMultiome_AllADVsCon_PFC_GabittoLabels/` |
| `MITMultiome_PFC_ROSMAP_AllADVsCon_GabLabels_runcode.R` | `…/MIT_ADMultiome_AllADVsCon_PFC_ROSMAP_GabittoLabels/` |

---

## 2. Old code (`older_projection_code/`)

Per-dataset / per-comparison drivers that the new `AllComparisons` scripts replace.
Grouped by subfolder; output = `save_dir1` module directory under
`…/greedy_march_pipeline_output/`.

### `batiuk/` — Batiuk 2022 SCZ
- `batiuk_processing.R`, `batiuk_COPA_runcode.R`, `batiuk_COPA_runcode_old.R` → `…/Batiuk_2022/`

### `brainSCOPE/` — bulk megaset & ROSMAP projections onto brainSCOPE
- `brainscope_COPA_runcode.R` → `…/brainSCOPE/`
- `brainscope_COPA_runcode_CMC.R` → `…/brainSCOPE_CMC/`
- `brainscope_COPA_runcode_MultiomeBrain.R` → `…/brainSCOPE_MultiomeBrain/`
- `brainscope_COPA_runcode_SZBD.R` → `…/brainSCOPE_SZBDMultiseq/`
- `rosmap_brainscope_COPA_runcode_CMC.R` → `…/brainSCOPE_CMC_ROSMAP/`
- `rosmap_brainscope_COPA_runcode_SZBD.R` → `…/brainSCOPE_SZBDMultiseq_ROSMAP/`

### `jorstad_lein/` — Lein DFC & MTG reference
- `LeinDFC_runcode.R`, `LeinDFC_runcode_new.R` → `…/LeinDFC/`
- `LeinDFC_runcode_corrected_CBscVI.R` → `…/leinDFC_CB_scVI/`
- `LeinDFC_runcode_uncorrected_PB.R` → `…/leinDFC_pb_uncorrected/`
- `LeinMTG_runcode.R`, `LeinMTG_runcode_new.R` → `…/LeinMTG/`

### `mathys/` — Mathys 2019/2023
- `Mathys_runcode.R` → `…/Mathys2019_ABIanno/`
- `mathys_2023.R` → preprocessing/annotation of MIT multiome h5ad (no module output; writes `_sif.csv`, commented)

### `MIT_Multiome/` — MIT AD Multiome, non-Gabitto labels
8 scripts: PFC/MTC × {AllADVsCon, earlyVsCon, lateVsEarly, APOE}
- → `…/MIT_ADMultiome_{AllADVsCon|earlyVsCon|earlyVLate|lateVsEarly|APOE}_{MTC|DFC|PFC}/`

### `morabito/` — Morabito
- `Morabito_runcode.R`, `Morabito_runcode_old.R` → `…/Morabito_ABIanno/` (old also `…/Morabito/`)

### `rosmapMods_COPA/` — ROSMAP-module projections (16 scripts)
con/AD modules projected onto Morabito, Mathys, SEAAD2024 (A9 & MTG, all comparisons),
MITMultiome (DFC & MTG), plus `rosmap_diffCoExpression.R`.
- → various `…/rosmap_AD/rosmap_AD_*` and `…/rosmap_diffCoExpress/*` dirs (one per script)

### `seaAD2024_COPA_A9/` — SEA-AD 2024 A9/DFC
Preprocessing + per-comparison COPA (AllADVsCon, earlyVsCon, lateVsEarly, apoe, MorabitoMatched)
- preprocessing → `…/SEAAD2024/`, splitting → `…/SEAAD2024_*/sn_summary_tables/`
- `seaAD2024_COPA_load_data.py` → writes `.mtx` next to source h5ad in datasets dir
- `seaAD2024_earlyVsLate_results_exploratory_analysis.R` → reads `euclidean_distances/*.qs`, writes exploratory PNGs to `…/SEAAD2024_earlyVsLate/`

### `seaAD2024_COPA_MTG/` — same as A9 but MTG region
- → `…/SEAAD2024_*_MTG/` dirs

### Top-level of `older_projection_code/`
- `COPA_conserve_seaAD2024_annotations.R` → `…/COPA_conserve/{mathys_vs_morabito, …_vs_SEAAD2024, morabito_vs_SEAAD2024}/`
- `dcopa_compare_runcode.R`, `dcopa_compare_pval_runcode.R` → `…/dcopa_shared/gabitto_vs_liu_DFC_*/`
- `yang2022_projection.R` → `…/Yang2022_ABIanno/`

---

## 3. Related analyses (`related_analyses/`)

Downstream / supporting analyses (Python projection-index calculators, plotting,
random-module null comparisons).

### `calc_projIndexSE_in_python/`
Python recompute of per-cell-type projection indices. Writes
`sn_proj_indices/.../indices_se_*.csv` into each dataset's module output dir.
- Active dataset scripts: `calc_projIndexSE_{brainSCOPE, brainSCOPE_v2…v7, lein_A9, lein_MTG, morabito}.py`, `..._mitMultiome_{MTC,PFC}[_*].py`, `..._mitMultiome_*_optimized_loop.py`, `..._SEAAD2024_{DFC,MTG}[_APOE/_ROSMAP]_optimized_loop.py`
- Drivers: `calc_projIndexSE_brainSCOPE.sh`, `run_all_projIndexSE.sh`
- `save_modlist_as_json.R` → writes `topmodposbc_table.json` / `seed_table.json` into `kme_tables/` of LeinDFC, rosmap, Brainseq_SCZ
- `old/` → `calc_projIndexSE_SEAAD2024_old.py`, `…_optimized.py` (superseded)

### `dotplots/` — summary dotplots / pairwise correlation figures
- `dotplots.R`, `dotplot_related_summaryPlots.R`, `dotplots_combined_pairwise_cors.R`, `fxns.R` → ggsave figures into `…/SEAAD2024_AllADVsCon_DFC/` and `/home/gugene/code/git/Consensus-analysis/Analyses/Projection_analyses/related_analyses/rand_mod_analysis/...`

### `random_mod_proj_analyses/` — real-vs-random module null comparisons
- `analysis_of_random_module_projections.R` → `/home/gugene/code/git/Consensus-analysis/Analyses/Projection_analyses/related_analyses/rand_mod_analysis/` (CSVs + plots)
- `realvsrand_dCoPA_analysis.R` → `…/rand_mod_dcopa_pvals/SEAAD2024/`
- `old_code/` → `realvsrand_dCoPA_analysis_{brainSCOPE,MIT}_old.R`, `realvsrand_summary_dotplots_old.R`, `realvsrand_summary_SEAAD2024_vs_MIT_old.R` (→ `figure_3/rand_mod_analysis/...`)

### Top-level of `related_analyses/`
- `calculate_rand_euclidean_distances.R` → `fwrite` random-distance tables
- `summary_of_COPA_results.R` → `~/…/copa_compare_summary_euler/` (Euler/summary plots)
- `umi_per_cell_3_datasets.R` → ggsave UMI/cell QC figure

---

## 4. Other / uncategorized

Within the folders there are **no orphaned data files** — every file is either a script
or tied to one. The only non-code artifacts checked into the tree:

- `calc_projIndexSE_in_python/logs/calc_projIndexSE_SEAAD2024_DFC_optimized_loop.log`
- `calc_projIndexSE_in_python/logs/calc_projIndexSE_SEAAD2024_MTG_optimized_loop.log`

Both are run logs of the two SEAAD2024 optimized-loop scripts.

**External dependency (not in this folder, but every script needs it):**
`/home/gugene/code/git/CoPA/` — `wrapper.R`,
`project_random_and_calculate_euclidean_standalone.R`,
`calculate_euclidean_distances_v2_standalone.R`.

---

## Cleanup flags

- **All real outputs live outside the repo** (`/mnt/bdata`, `~/figure_sup_analyses`), so this
  folder is code-only and not self-contained without those paths.
- The new `*_AllComparisons*` scripts make most of `older_projection_code/MIT_Multiome/`,
  `seaAD2024_COPA_*/`, and the `deprecated/` set redundant — that is the clearest dead-weight
  if cleaning up.

# SCZ gene database (fig_8)

Curated schizophrenia-linked gene database for figure 8 (SCZ↔AD shared modules),
DFC brain region. These CSVs are **tracked** (force-added past the `figures/**/*.csv`
gitignore) because they are the product of a **non-deterministic, agent-driven
literature-review pipeline** (LLM agents + PubMed/OMIM/Open Targets/ClinVar
lookups) and are not cheaply or reproducibly regenerable — they are effectively
primary curated inputs the figure depends on.

## Pipeline (produced by `scz_db.R` + the `*.py` helpers here)

1. Extract unique candidate genes for DFC (119 pipeline-track genes from `dfc_overlaps.csv`).
2. Prefilter against OMIM (MIM 181500), Open Targets (MONDO:0005090, curated
   sources only), ClinVar — `build_scz_prefilter.py`.
3. Three-agent review: Agent 1 PubMed search → Agent 2 PMID verification →
   Agent 3 reconciliation.
4. Merge pipeline + prefilter tracks into the summary table.
5. Build `scz_db.R` (`scz_db_dfc` R list) + consistency checks.

## Files

| File | Contents |
|------|----------|
| `scz_gene_databases_prefilter.csv` | Prefilter DB: 2,196 SCZ-associated genes (OMIM 268, OpenTargets 1,850, ClinVar 285). |
| `scz_db_agent1_dfc.csv` (+ `_batch_1..5`) | Agent 1 PubMed-search output: 22 pipeline genes retained (of 119; 44 had ≥1 hit). Batches are the chunked agent runs. |
| `scz_db_agent2_dfc.csv` | Agent 2 PMID verification (all 22 verified; 10 "Y with note"). |
| `scz_db_agent3_dfc_pipeline.csv` | Agent 3 reconciliation (all 22 confirmed). |
| `scz_db_summary_table_dfc.csv` | Collated: **56 genes** (22 pipeline + 34 prefilter, no overlap). Read by `fig8.R`. |
| `scz_genes.txt` | Reference gene list (10 known SCZ genes) used for cross-checks. |
| `scz_db.R`, `scz_ref_index.R` | Code that builds/uses `scz_db_dfc` (already tracked). |

Final `scz_db_dfc` = 56 entries. Pipeline genes carry integer PMID refs +
paper counts; prefilter genes carry `ref="website"` (OMIM 2, ClinVar 6,
OpenTargets 26). Two false-positive symbols excluded (GK, OAT — matched clinical
abbreviations, not genes).

*Pipeline completed 2026-04-29.*

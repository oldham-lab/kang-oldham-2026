# fig_7 — AD/SCZ gene-database agent artifacts

Primary records of the LLM-assisted AD and SCZ gene curation behind Table S15
and the fig_7 gene-database panels. Consolidated here from the superseded
fig_7/v7.2 and fig_8/v1 folders so the provenance sits beside the current
build; the artifacts themselves are unchanged. These are **not regenerable**: the pipeline is
non-deterministic, so re-running it produces a different result. They are
tracked past the `figures/**/*.txt` ignore for that reason.

| File | What it is |
|---|---|
| `prompt_ad_db_v2.txt` | Prompt driving the v2 run |
| `prompt_ad_db_v3.txt` | Revised prompt driving the v3 MTG run |
| `prompt_changes.txt` | Summary of what changed between v2 and v3 |
| `alzheimers_genes_ranked.txt` | 19 well-established AD risk genes, 4 tiers |
| `mtg_v3_candidates.txt` | 403 MTG candidate genes for the v3 run |
| `ad_gene_databases_citations.txt` | Citations for the source gene databases |
| `opentargets_datasources.txt` | Descriptions of all 20 OpenTargets datasources |
| `ad_db_summary_table_dfc_opentargets_sources.txt` | Per-gene OpenTargets evidence, DFC |
| `ad_db_summary_table_mtg_opentargets_sources.txt` | Per-gene OpenTargets evidence, MTG |
| `prompt_scz_db_v1.txt` | Prompt driving the three-agent SCZ run (superseded by fig_8/v4) |

Methods narrative: `../ad_scz_gene_db_methods.md`. The current SCZ pipeline and
its own artifacts live in `../../fig_8/v4/scz_db_pipeline/`.
The superseding v8 pipeline lives in `../../v8/ad_db_pipeline/`.

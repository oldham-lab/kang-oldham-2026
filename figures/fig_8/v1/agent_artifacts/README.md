# fig_8 v1 — SCZ gene-database agent artifacts

Primary records of the first (three-agent) LLM-assisted SCZ gene curation.
**Not regenerable** — the pipeline is non-deterministic — so these are tracked
past the `figures/**/*.txt` ignore.

| File | What it is |
|---|---|
| `prompt_scz_db_v1.txt` | Prompt driving the v1 three-agent SCZ run |
| `dfc_pipeline_track.txt` | Genes on the literature-pipeline track (DFC); read by `../fetch_pubmed_scz.py` |
| `dfc_prefilter_track.txt` | Genes on the prefilter track (DFC) |

The curated CSV outputs of this run (`../scz_db_*.csv`,
`../scz_gene_databases_prefilter.csv`) are tracked in the parent directory and
documented in `../README_scz_db.md`. The superseding v4 pipeline lives in
`../../v4/scz_db_pipeline/`.

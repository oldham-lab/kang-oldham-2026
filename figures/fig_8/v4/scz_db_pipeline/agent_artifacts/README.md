# fig_8 v4 — SCZ gene-database agent artifacts

Primary records of the HGNC-corrected SCZ gene-database rebuild that supersedes
the fig_8 v1 three-agent run. **Not regenerable** — the classification step is
non-deterministic — so these are tracked past the `figures/**/*.txt` ignore.

| File | What it is | Read by |
|---|---|---|
| `prompt_scz_classify.txt` | Classification prompt for the rebuild | — |
| `dfc_candidates.txt` | Candidate genes considered (DFC) | `../assemble_scz_db.py` |
| `literature_track.txt` | Genes on the literature track | `../alias_fetch.py`, `../ncbi_prefetch.py` |
| `prefilter_track.txt` | Genes on the prefilter track | `../prefilter_counts.py` |

The pipeline's JSON inputs, caches and classification results sit in the parent
directory and are tracked already. Methods narrative:
`../../../fig_7/v8/ad_scz_gene_db_methods.md`.

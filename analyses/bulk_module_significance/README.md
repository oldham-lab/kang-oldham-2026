# Bulk module significance filter

Small per-module significance counts derived from the bulk pairwise-correlation
analysis, used as a **module pre-filter** across many figures (fig_3, fig_4,
fig_6, fig_7, fig_7_sup) and the Shiny precalc: a module is kept when it has at
least 2 significant bulk correlations (`vals >= 2`).

| File | Meaning |
|------|---------|
| `bulk_cors_sigcount_bonf_1158.csv` | Bonferroni significant-correlation counts per module (1158-module set). |
| `bulk_cors_sigcount_bonf_rosmapAD_1127.csv` | Same, ROSMAP-AD variant (1127-module set). |

These were originally produced under `Code_for_figures/fig_3_old/bulk_cor_significance_analysis/`
(now in the `consensus-analysis-archive` repo) and were previously untracked.
They are relocated here and tracked because they are live reproduction inputs,
not generated outputs of the current figure scripts.

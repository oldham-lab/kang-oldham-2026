# fig_1 DE-overlap gene-symbol vintage bug — HANDOFF (logged 2026-07-13)

Self-contained note for a future session (possibly a different account) to fix.
**Status: CODE FIX APPLIED 2026-07-13 — DE/figure REGENERATION STILL PENDING.**
The canonicalization is implemented in `full_DE_pipeline_optimized_v2.4.R`
(load-time; see "Fix as implemented" below). Still TODO: re-run the DE pipeline
(STEP 3 onward) to regenerate the `DE_v2.4/{DESeq2,edgeR}` `.qs` objects, then
re-run `v3/fig1_v3.R` and report the overlap delta.

Discovered while auditing figure impact of the
HGNC harmonization fix applied to the *projection* pipelines. This is a
**separate, pre-existing** bug in fig_1's DE pipeline — the projection fix
(COPA repo, commit `b94a6b3`) did NOT touch it.

## Fix as implemented (2026-07-13)

Load-time canonicalization (does NOT require re-pseudobulking — the STEP 1/2
CSVs are unchanged). Added to `full_DE_pipeline_optimized_v2.4.R`:
- `ALIAS_MAP_PATH` config + `canonicalize_symbols()` / `read_pb_canonical()`
  helpers (the latter renames gene symbols to current HGNC and **sum-collapses**
  post-rename collisions so rownames stay unique).
- `common_genes_by_region` now intersects canonicalized symbols.
- STEP 3/4 (edgeR) and STEP 5 (DESeq2) load pseudobulk via `read_pb_canonical`.
Verified: symbols unique, sample cols preserved, total counts conserved, row
drop = collision count (Jorstad DFC −30, SEAAD DFC −110), universe +226/+231.
Because the pseudobulk CSVs are untouched, the regeneration can start at the
`common_genes_by_region` block (STEP 3) and skip the expensive STEP 1/2.

## The bug (one paragraph)

fig_1's differential-expression gene universe is the **symbol-string
intersection** of Jorstad (Lein) and SEAAD (Gabitto) pseudobulk matrices. But
the two datasets use different HGNC symbol vintages — **Jorstad = current, SEAAD
= older** — so genes renamed between vintages (AARS↔AARS1, HIST1H2BJ↔H2BC11,
MARCH1↔MARCHF1, histones, aminoacyl-tRNA synthetases, etc.) fail to match and
are **dropped from DE testing in BOTH datasets** before DESeq2/edgeR run. They
therefore cannot appear in either dataset's DE set or in the Jorstad∩SEAAD
"shared DE genes" overlap that fig_1 reports — the overlap is undercounted.

## Location

`Code_for_figures/fig_1/full_DE_pipeline_optimized_v2.4.R`
- **lines ~270–276**: `common_genes_by_region <- ... intersect(jg, sg)` where
  `jg`/`sg` are column-1 gene symbols of
  `{jorstad,seaad}_cell_expression_by_donor_subclass_sum_<REGION>.csv`
  (staged under `DE_v2.4/pseudobulk/`). This defines the per-region DE universe.
- **line ~391**: a second `intersect(rownames(jorstad), rownames(sea))` in the
  combined-DE helper — same fix needed.

Overlap consumed in `v3/fig1_v3.R`: euler (`~174–177`), shared-DE intersection
(`~292–294`: `intersect(edgeR_jorstad, edgeR_sea)` ∩ `intersect(DESeq2_jorstad,
DESeq2_sea)`), and the UpSet. `recode_map` in fig1_v3.R is a **cell-type** rename
only — no gene-symbol handling anywhere.

## Quantified impact (measured 2026-07-13 vs the alias map)

Genes recoverable into the DE universe if both sides are canonicalized to
current HGNC before intersecting:

| region | current DE universe | harmonized | recovered |
|---|---|---|---|
| DFC | 20,342 | 20,568 | **+226** |
| MTG | 20,661 | 20,892 | **+231** |

+226/+231 is the CEILING on genes entering DE; the actual change to the *shared*
overlap is the subset of those that are DE-significant in both datasets
(unknowable without re-running DE).

## Fix approach

Alias map (prev_symbol → current_symbol, 56,964 entries):
`/home/gugene/code/git/COPA/dataset_processing_unified/hgnc_prev_to_current.tsv`
(built by `build_alias_map.R`; from HGNC complete set).

1. Before the intersect, **canonicalize the rownames** of both pseudobulk
   matrices to current HGNC (map `prev`→`current`; leave unmapped as-is).
   This RENAMES rows (not just adds match keys, unlike the projection fix)
   because DE needs the same symbol for the same gene in both datasets and
   DESeq2/edgeR require **unique** rownames.
2. Handle collisions: if canonicalization makes two rows share a symbol (or a
   current symbol already present), decide a dedup rule (sum counts, or
   first-wins) — check how many collisions occur first; expected to be few.
3. Intersect the canonicalized symbols → new common universe.
4. Re-run `full_DE_pipeline_optimized_v2.4.R` to regenerate
   `DE_v2.4/{DESeq2,edgeR}/{jorstad,seaad}/*_matchedGenes.qs` (DFC + MTG).
5. Re-run `v3/fig1_v3.R` (+ `v3/supplementary.R`) to rebuild the figure; report
   the overlap delta vs the published version.

## Notes / gotchas

- fig_1's other inputs (bulk fidelity, gene means, cell annotations) are
  unaffected — only the Jorstad↔SEAAD DE intersection matters.
- Related but distinct: the projection-side fix is in the COPA repo
  (`HGNC_HARMONIZATION_HANDOFF.md`, commit `b94a6b3`); Morabito + the three
  Python datasets were regenerated there. This fig_1 DE issue is in the
  Consensus-analysis repo and stands alone.
- Watch fig_1 object-ordering gotcha (DE objects stored in non-cts-subclass
  order) when re-running — see prior fig_1 notes.

# fig_7 incremental AI literature review — scope (run IN v8/, minimal-file output)

> **⚠️ SUPERSEDED PLANNING NOTE.** This describes an early *incremental* plan that was NOT
> what shipped. The final fig_7 build was a **FULL rerun** on all overlap genes, and
> **DFC used BOTH tracks** (literature + curated-DB prefilter), not literature-only.
> Actual output: `ad_db_summary_table_dfc.csv` = 74 genes (39 literature + 3 curator +
> **32 prefilter**). See `ad_scz_gene_db_methods.md` for the authoritative as-built method (now covering both AD and SCZ); trust
> the summary-table output over this note.

Per prompt_ad_db_v3.txt. DFC = literature-pipeline ONLY (no prefilter track). MTG = both tracks.
Only genes NEWLY gained in the post-HGNC overlap (not in the existing ad_db_summary_table) are reviewed.

## Literature-pipeline genes (Agent 1 search -> Agent 2 verify PMIDs -> Agent 3 reconcile)
DFC (all 19 gained): ATP2A2, CNOT2, CTDSPL2, ECPAS, KIAA0586, KLHL20, LRPPRC, MTERF2, NAA35,
  NRDC, PGAP4, SEZ6L, TMEM260, TSGA10, APLP2, MTMR2, SCAPER, UBE3B, UBE4B
MTG (8 novel, not in prefilter): ATP6V1B2, CNOT2, MARCHF7, NAA35, NECAP1, PEX19, SLC12A6, TMEM260
Unique genes to search (CNOT2/NAA35/TMEM260 shared): 24
  -> for shared genes use identical n/ref/mechanism in both ad_db_dfc and ad_db_mtg (prompt Step 4).

## Prefilter-track (mechanical add, MTG only, no lit search)
COG5  (in ad_gene_databases_prefilter.csv)

## Method (prompt Steps 3,4,4.5,4.6)
- Agent1: NCBI esearch/efetch, query `"GENE"[TIAB] AND "Alzheimer"[TIAB]`, 2000-present, top 5, ~0.4s delay.
  Include only genes with a verified PMID linking to one of the 13 mechanism categories. False-positive symbol check.
- Agent2: independently verify each PMID (real + on-topic); substitute/flag.
- Agent3: reconcile; include if >=2 agents agree evidence is supported.
- Step 4.5: Pubmed_total_hits (esearch Count) for each new gene.
- Step 4.6: OpenTargets AD score (MONDO_0004975) for each new gene (GraphQL; correct query shape TBD).

## Minimal-file OUTPUT (touch only these; no agent-batch/checkpoint sprawl)
1. append passing genes (+COG5) as rows to v8/ad_db_summary_table_dfc.csv and v8/ad_db_summary_table_mtg.csv
   (match existing schema: Gene, Region, N_papers, AD_mechanism, References, References_verified,
    AD_involvement_verified, Agent3_resolution, Pubmed_total_hits, OpenTargets_AD_score)
2. regenerate v8/ad_db.R (ad_db_dfc, ad_db_mtg)
Then SYNC these 3 files v8/ -> canonical v7.2/ (fig_6 + table_s15 read from v7.2).

## APIs verified reachable: NCBI E-utilities (APLP2=137 hits), OpenTargets GraphQL (needs correct field path).

## After review: rerun v8 pipeline (regen already done -> fig7_v7.1.R -> fig7_v7.2_summary_table.R -> assemble_figure.py),
## verify figure, check fig_6 ripple (fig6_v4.R reads ad_db_summary_table).

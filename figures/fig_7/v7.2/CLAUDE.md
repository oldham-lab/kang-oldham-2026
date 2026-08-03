# AD Gene Database Pipeline — v7.2 Context

## What this directory is

Output of a multi-agent pipeline that identifies Alzheimer's disease (AD)-linked genes
from brain region overlap gene lists, verifies them against the published literature,
and produces an R list object (`ad_db.R`) for use in figure 7 code.

Two brain regions are processed: **DFC** (dorsolateral prefrontal cortex) and **MTG**
(middle temporal gyrus).

---

## Input files

- `dfc_overlaps.csv` — 198 unique genes from DFC overlap analysis
- `mtg_overlaps.csv` — 544 unique genes from MTG overlap analysis
- `prompt_ad_db_v2.txt` — the prompt that drove the original pipeline
- `prompt_ad_db_v3.txt` — revised prompt driving the v3 MTG run (see below)

---

## Pipeline structure (from prompt_ad_db_v2.txt)

**Step 1** — Extract unique genes per region from overlap CSVs.

**Step 2** — Pre-filter against curated AD gene databases (DisGeNET C0002395, OMIM,
AlzGene) using bulk queries. Only genes with an entry in at least one database proceed.

**Step 3** — Three-agent literature review pipeline (DFC first, then MTG):
- Agent 1: PubMed literature search; only outputs genes with a direct AD mechanism
  link (12 categories) supported by at least one verified PMID.
- Agent 2: Verifies up to 3 PMIDs per gene; flags discrepancies with qualified-Y notes.
- Agent 3: Reconciles flagged genes using majority-vote inclusion rules.

**Step 4** — Build `ad_db.R` with `ad_db_dfc` and `ad_db_mtg` R list objects.

**Step 5** — Consistency checks across all outputs.

---

## v3 prompt changes (prompt_ad_db_v3.txt)

Created 2026-04-27. Key differences from v2:

1. **Inverted prefilter logic** — Step 3 now runs on genes NOT in the prefilter
   database (novel candidates), rather than genes that ARE in it.
2. **Single-gene PubMed search** — multi-gene OR queries removed; one gene at a time,
   top 5 results sorted by relevance (was top 10).
3. **Agent 2 verifies all 5 PMIDs** — was top 3 most-cited only.
4. **13th mechanism category added** — "Implication in a genome-wide association
   study (GWAS) of AD."
5. **N_papers capped at 5** — "5+" overflow notation (was "10+").

---

## v3 MTG run — status (as of 2026-04-28)

### Scope
- 403 candidate genes: MTG overlap genes not in prefilter DB and not already in
  `ad_db_summary_table_mtg.csv`.
- Candidate list saved to `mtg_v3_candidates.txt` (sorted alphabetically).
- 17 batches of 25 (final batch: 5 genes).

### Agent 1 progress — COMPLETE
| Batch | Genes | Reviewed | Included | Status |
|---|---|---|---|---|
| 1 | AASDHPPT–ARL15 | 25 | 5 | Complete |
| 2 | ARL6IP5–CADM3 | 25 | 1 | Complete |
| 3 | CAMSAP1–CSTF2 | 25 | 6 | Complete |
| 4 | CUL1–ELOVL6 | 25 | 8 | Complete |
| 5 | EMC1–GPALPP1 | 25 | 7 | Complete |
| 6 | GPALPP1–KLHL7 | 25 | 5 | Complete |
| 7 | KLHL9–MAP1LC3B | 25 | 2 | Complete |
| 8 | MAP2K1–NDFIP1 | 25 | 4 | Complete |
| 9 | NDFIP2–PIK3CA | 25 | 8 | Complete |
| 10 | PIK3C2A–RANBP9 | 25 | 8 | Complete |
| 11 | RAP1GDS1–RORB | 25 | 1 | Complete |
| 12 | RTN1–SKP1 | 25 | 8 | Complete |
| 13 | SLC9A6–TAOK1 | 25 | 5 | Complete |
| 14 | TFB2M–TSPAN5 | 25 | 3 | Complete |
| 15 | TTC39C–VPS41 | 25 | 7 | Complete |
| 16 | VPS50–ZNF654 | 25 | 4 | Complete |
| 17 | ZNF670–ZXDA | 5 | 0 | Complete |
| **Total** | | **403** | **82** | **Complete** |

Merged output: `ad_db_agent1_mtg_v3.csv` (82 genes, 83 rows including header).

### Agent 2 — COMPLETE (2026-04-28)
Output: `ad_db_agent2_mtg_v3.csv` (82 genes). All genes verified Y/Y.
6 qualified flags: OXR1 (corrected, false-positive PMID removed), SEC31A, HECW2, PPP2R2C, OPA1, PRKAA2.

### Agent 3 — COMPLETE (2026-04-28)
Output: `ad_db_agent3_mtg_v3.csv` (82 genes). All 82 confirmed.
Qualified genes resolved with justification notes in Agent3_resolution column.

### Integration — COMPLETE (2026-04-28)
- `ad_db_summary_table_mtg.csv`: 222 genes (126 pipeline-only + 18 pipeline/prefilter overlap + 78 prefilter-only)
- `ad_db_summary_table_mtg_pipeline_only.csv`: backup of 144-gene pipeline-only table (CLDN11 removed)
- `ad_db.R`: ad_db_mtg has 222 entries; prefilter-only genes use n=1 and string ref (website); 18 overlap genes have refs updated to database websites
- `ad_implicated_genes.txt`: deleted
- CLDN11: removed from all MTG outputs — was in prefilter DB and included by v3 Agent 1, but is not in mtg_overlaps.csv

### Search methodology notes
- PubMed E-utilities API (esearch + efetch); query: `"GENE"[TIAB] AND "Alzheimer"[TIAB]`
- Top 5 results by relevance; 2000–present; abstracts only (not full text).
- Rate: ~0.4 s delay between queries to respect NCBI rate limits.

---

## Output files and what they contain

| File | Contents |
|---|---|
| `mtg_v3_candidates.txt` | 403 MTG candidate genes for v3 run (sorted alphabetically) |
| `ad_db_agent1_mtg_v3_batch_{1-17}.csv` | v3 Agent 1 batch checkpoints (all 17 batches complete) |
| `ad_db_agent1_mtg_v3.csv` | v3 Agent 1 merged output: 82 novel AD-linked MTG genes |
| `ad_db_agent2_mtg_v3.csv` | v3 Agent 2 verification: 82 genes, all verified Y (6 qualified) |
| `ad_db_agent3_mtg_v3.csv` | v3 Agent 3 reconciliation: 82 genes confirmed |
| `ad_db_agent1_dfc.csv` | Agent 1 output: 56 AD-linked DFC genes |
| `ad_db_agent1_dfc_batch_{1-3}.csv` | Output batches (25/25/6 rows); exactly match merged file |
| `ad_db_agent2_dfc.csv` | Agent 2 verification of DFC genes |
| `ad_db_summary_table_dfc.csv` | Final DFC table (77 genes: 50 pipeline-only + 6 pipeline/prefilter overlap + 21 prefilter-only) |
| `ad_db_summary_table_dfc_pipeline_only.csv` | Backup of 56-gene pipeline-only DFC table (PMID refs intact for the 6 overlap genes) |
| `ad_db_agent1_mtg.csv` | Agent 1 output: 63 AD-linked MTG genes (merged from two runs) |
| `ad_db_agent1_mtg_batch_{1-3}.csv` | Batch files from the supplementary new MTG run only (20 genes) |
| `ad_db_agent1_mtg_new.csv` | Supplementary MTG run Agent 1 output (20 genes, all in merged file) |
| `ad_db_agent2_mtg.csv` | Agent 2 verification of full MTG gene set (63 genes) |
| `ad_db_agent2_mtg_new.csv` | Agent 2 output for supplementary run (20 genes) |
| `ad_db_agent3_mtg_new.csv` | Agent 3 output for supplementary run (20 genes) |
| `ad_db_summary_table_mtg.csv` | Final MTG table (222 genes: 126 pipeline-only + 18 pipeline/prefilter overlap + 78 prefilter-only) |
| `ad_db_summary_table_mtg_pipeline_only.csv` | Backup of 144-gene pipeline-only MTG table (PMID refs intact for the 18 overlap genes; CLDN11 removed) |
| `ad_db.R` | Final R list objects: `ad_db_dfc` (77 genes) and `ad_db_mtg` (222 genes) |
| `prompt_changes.txt` | Notes on prompt iterations |

---

## Known pipeline execution history

### Step 2 (prefilter) was skipped in all runs
The prompt specifies a database prefilter before Agent 1. This was not applied.
Evidence: 50 DFC and 45 MTG final genes have no entry in OMIM, Open Targets, ClinVar,
or AlzForum — they could not have passed a correctly applied prefilter.
Agent 1 acted as the de facto filter, reviewing the full overlap gene lists and
only outputting genes where it found a PubMed AD link.

### DFC run history
Earlier run(s) were prematurely terminated (no surviving batch files). A final
complete run processed all 198 DFC overlap genes and produced 56 AD-linked genes,
written as output batches of 25/25/6. The DFC pipeline is complete.

### MTG run history
Earlier run(s) were prematurely terminated. Surviving output from those runs
contributed 43 genes to `ad_db_agent1_mtg.csv` (batch files not preserved).
A supplementary new run was then launched, adding 20 more genes
(`ad_db_agent1_mtg_new.csv`, current batch files). Combined total: 63 genes.
The MTG pipeline is complete but covers only 63 of 544 overlap genes.

---

## Verification findings (from session analysis, 2026-04-27)

### Agent flagging
- No genes received a hard "N" flag from Agent 2 or 3.
- DFC: 14 genes received qualified-Y flags (indirect/pathway-level AD links).
- MTG: 4 genes received qualified-Y flags (preprint PMIDs; JAZF1 genetic-only link).
- All flagged genes were confirmed by Agent 3.

### Final table completeness
- `ad_db_summary_table_mtg.csv` contains all 63 genes from all MTG runs. ✓
- All genes in both summary tables are present in `ad_db.R`. ✓

---

## Prefilter database (ad_gene_databases_prefilter.csv)

A reference table of 2,257 AD-associated genes scraped from publicly accessible
curated databases (generated 2026-04-27). Columns: Gene, Gene_name, In_OMIM,
In_OpenTargets_Platform, OpenTargets_score, OpenTargets_sources, In_ClinVar,
In_AlzForum, Databases.

Sources successfully queried:
- **OMIM** — 362 genes (via NCBI E-utilities)
- **Open Targets Platform** (includes merged Open Targets Genetics) — 1,913 genes
  (MONDO_0004975; curated datasources only)
- **ClinVar** — 97 genes (via NCBI E-utilities; 1,177 AD variants)
- **AlzForum** — 6 genes (partial; 403 blocked full scrape)

Sources inaccessible: DisGeNET (API key required), AlzGene (defunct), GWAS Catalog
(download endpoints returning 404), AGORA (API not accessible), MalaCards/GeneCards
(403 Forbidden), Ensembl BioMart (no results for AD queries).

See `ad_gene_databases_citations.txt` for full APA citations.

### Prefilter coverage of overlap genes
| Region | Overlap genes | In prefilter DB | Absent from DB |
|---|---|---|---|
| DFC | 198 | 27 (14%) | 171 (86%) |
| MTG | 544 | 96 (18%) | 448 (82%) |

### Genes that passed prefilter but were missed by pipeline
- **DFC: 21 genes** — reviewed by Agent 1 but found no qualifying AD evidence.
  Notable: TUBA1B, TUBA4A, TUBB, TUBB2A, TUBB3, TUBB4B (tubulin family), VPS13C.
- **MTG: 78 genes** — most never reviewed due to premature termination of earlier runs.
  Notable missed genes with strong AD relevance: C9orf72, MME (neprilysin), VPS35,
  TOLLIP, TBK1, CACNA1B, and the full tubulin family.

---

## OpenTargets AD association score annotation (2026-05-21)

All genes in both summary tables were annotated with their Open Targets Platform
AD association score (disease: MONDO_0004975, data version 26.03) using the
GraphQL API (https://api.platform.opentargets.org/api/v4/graphql).

### Method
- Gene symbols resolved to Ensembl IDs via the `search` query
- AD association scores retrieved via `associatedDiseases(BFilter: "MONDO_0004975")`
  in batches of 20 using GraphQL aliases
- Datasource-level evidence retrieved via `evidences(efoIds: ["MONDO_0004975"])`
  for all scored genes; aggregated by datasource (max score per datasource per gene)
- Script: /tmp/ot_query_genes.py (MTG) and /tmp/ot_query_dfc.py (DFC)

### Results

| Region | Total genes | With score | No score |
|--------|-------------|------------|----------|
| MTG    | 222         | 214        | 8        |
| DFC    | 77          | 74         | 3        |

Genes with no AD score in OpenTargets:
- MTG: EDEM3, FASTKD2, GSE1, RAB11FIP5, RBM11, RNF11, RPL39, SEC24A
- DFC: EDEM3, RPL32, RPL39

Score range (scored genes): ~0.001 (ribosomal proteins) to ~0.52 (NDUFS1 in MTG).
Top MTG scorers: NDUFS1 (0.523), PIKFYVE (0.509), MTOR (0.494).

### Files modified / created

| File | Change |
|------|--------|
| `ad_db_summary_table_mtg.csv` | `OpenTargets_AD_score` column added (last column) |
| `ad_db_summary_table_dfc.csv` | `OpenTargets_AD_score` column added (last column) |
| `ad_db_summary_table_mtg_opentargets_sources.txt` | New — datasource evidence per scored gene |
| `ad_db_summary_table_dfc_opentargets_sources.txt` | New — datasource evidence per scored gene |
| `opentargets_datasources.txt` | New — general descriptions of all 20 OpenTargets datasources |
| `prompt_ad_db_v3.txt` | Step 4.6 added documenting the OpenTargets annotation procedure |

### Notes on interpretation
- The overall score is a harmonic sum across datasource scores, capped at 1.
- Most pipeline genes are implicated via `europepmc` (literature text mining) only —
  this reflects co-mention frequency, not curated biological evidence.
- Genes with `gwas_credible_sets` as a datasource have the strongest genetic evidence.
- The ribosomal protein genes (RPL/RPS family) score very low (0.001–0.01),
  consistent with indirect AD evidence in the literature.
- NRXN1's top evidence is a GWAS credible set from a PRS-derived phenotype study
  (GCST90132260, Gouveia 2022), which is methodologically indirect.
- VPS35's AD link is almost entirely literature co-mention; its primary association
  in OpenTargets is with Parkinson disease.

---

## Cross-check: ranked classic AD genes vs pipeline outputs (2026-04-29)

`alzheimers_genes_ranked.txt` contains 19 well-established AD risk genes across 4 tiers
(APOE, PSEN1, APP, PSEN2, TREM2, SORL1, ABCA7, BIN1, CLU, PICALM, CR1, CD33, MS4A4A,
MS4A6A, EPHA1, INPP5D, MEF2C, FERMT2, HLA-DRB1).

All 19 are present in `ad_gene_databases_prefilter.csv`. Subset reaching final outputs:

| Gene | ad_db_summary_table_dfc | ad_db_summary_table_mtg | ad_db.R |
|---|---|---|---|
| APP | yes | yes | yes (both) |
| PSEN1 | — | yes | yes (MTG) |
| TREM2 | — | yes | yes (MTG) |
| SORL1 | — | yes | — |
| ABCA7 | — | yes | — |
| BIN1 | — | yes | — |

Remaining 13 (APOE, PSEN2, CLU, PICALM, CR1, CD33, MS4A4A, MS4A6A, EPHA1, INPP5D,
MEF2C, FERMT2, HLA-DRB1) appear only in the prefilter database — they were not in
the DFC or MTG overlap gene lists from the brain-region analyses.

---

## SCZ Gene Database Pipeline — v1 (fig_8/v1)

A parallel schizophrenia pipeline was run for figure 8, covering the DFC brain region.
Located at `../../../fig_8/v1/` relative to this directory.

### Overview

Produces `scz_db.R` (`scz_db_dfc` R list object) from DFC overlap gene list.
Mirrors the AD pipeline structure with the same three-agent review approach.

### Input

- `dfc_overlaps.csv` — 119 pipeline-track candidate genes
- `prompt_scz_db_v1.txt` — pipeline prompt

### Prefilter database (`scz_gene_databases_prefilter.csv`)

2,196 SCZ-associated genes from curated sources (built by `build_scz_prefilter.py`):
- **OMIM** — 268 genes (MIM 181500 via NCBI E-utilities)
- **Open Targets Platform** — 1,850 genes (MONDO:0005090, curated datasources only)
- **ClinVar** — 285 genes (≤10 genes per variant; filters out large CNVs)
- **SZDB** — 0 genes (www.szdb.org inaccessible — connection reset)

### Pipeline run — COMPLETE (2026-04-29)

**Agent 1** (`scz_db_agent1_dfc.csv`, 5 batch files):
- 119 candidates reviewed; 44 had ≥1 PubMed hit; 22 included
- Query: `"GENE"[TIAB] AND "schizophrenia"[TIAB]`, top 5 by relevance, 2000–present
- Two false-positive symbols excluded: GK (matched clinical abbreviation),
  OAT (matched "oral anticoagulant therapy")

**Agent 2** (`scz_db_agent2_dfc.csv`):
- All 22 genes verified; 10 received "Y with note" qualifications
- All 60 PMIDs confirmed real via `pmid_verification.json`

**Agent 3** (`scz_db_agent3_dfc_pipeline.csv`):
- All 22 confirmed (rule 1 majority-vote for all qualified cases)
- No hard "N" flags from any agent

**Integration** (`scz_db_summary_table_dfc.csv`):
- 56 genes total: 22 pipeline + 34 prefilter; no overlap between tracks
- Prefilter breakdown: OMIM 2, ClinVar 6, OpenTargets_Platform 26

**Final output** (`scz_db.R`, `scz_db_dfc`):
- 56 entries; pipeline genes use integer PMID refs and n = paper count ('1'–'5+')
- Prefilter genes use n=1 and string ref ("omim.org" / "platform.opentargets.org" /
  "ncbi.nlm.nih.gov/clinvar"), same convention as AD pipeline

### Cross-check: scz_genes.txt (10 canonical SCZ genes)

File: `scz_genes.txt` — SETD1A, GRIN2A, TRIO, CACNA1G, SP4, RB1CC1, CUL1, XPO7,
HERC1, SRRM2.

- **RB1CC1** — only gene present in both `dfc_overlaps.csv` and `scz_genes.txt`;
  entered final output via the prefilter track (not the literature pipeline)
- **SETD1A, GRIN2A, CACNA1G, SP4** — in prefilter DB but not in `dfc_overlaps.csv`
- **TRIO, CUL1, XPO7, HERC1, SRRM2** — not found in any pipeline file

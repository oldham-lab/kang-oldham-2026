# Disease gene-database construction — methods (Fig. 7 AD & Fig. 8 SCZ)

Rebuilt 2026-07-22. Describes how the disease-gene databases that annotate the reproducible
dCoPA gene lists were built. **One pipeline** produces both:

- **Alzheimer's disease (AD)** — Fig. 7 panel G / Table S15: `fig_7/v8/ad_db.R`,
  `ad_db_summary_table_{dfc,mtg}.csv`.
- **Schizophrenia (SCZ)** — Fig. 8 panel E / Table S15: `fig_8/v4/scz_db.R`,
  `scz_db_summary_table_dfc.csv`.

The two runs share the identical procedure and differ only in the disease query term, the
mechanism-category set, the curated prefilter database, and the input gene lists. This
supersedes the earlier three-agent procedure (`prompt_ad_db_v3.txt` for AD,
`prompt_scz_db_v1.txt` for SCZ).

**Disease-specific parameters**

| | AD (Fig. 7) | SCZ (Fig. 8) |
|---|---|---|
| Regions | DFC + MTG | DFC |
| dCoPA consensus | Gabitto + Liu single-nucleus | brainSCOPE CMC + SZBDMulti-Seq |
| Input overlap genes | 215 (DFC), 542 (MTG) | 151 (DFC) |
| Disease PubMed query | `"SYM"[TIAB] AND "Alzheimer"[TIAB]` | `"SYM"[TIAB] AND "Schizophrenia"[TIAB]` |
| Curated prefilter DB | OMIM, Open Targets, ClinVar, AlzForum (2,257 genes) | OMIM, Open Targets, ClinVar (2,196 genes; SZDB inaccessible) |
| Mechanism categories | 13 (AD) | 14 (SCZ) |
| Open Targets score | yes (MONDO_0004975) | not computed |
| **Final database** | **DFC 74, MTG 215** | **DFC 51** |

---

## Detailed methods (shared pipeline)

**Input.** Each database annotates the reproducible cross-dataset dCoPA gene lists — genes
belonging to modules significant in **both** single-nucleus projections of the disease's two
datasets (AD: Gabitto **and** Liu; SCZ: brainSCOPE **CMC** and **SZBDMulti-Seq**), separately
per region. AD: `dfc_overlaps.csv` / `mtg_overlaps.csv` = 215 / 542 unique genes. SCZ:
`dfc_overlaps.csv` = 151 unique genes. All gene symbols are HGNC-harmonized before the
overlap. Each gene enters via one of three tracks.

**1. Track assignment (both regions, both diseases).** Genes present in the disease's curated
prefilter database are added directly (**prefilter track**), with the database as their
reference and no literature search; the remaining genes go through the literature pipeline.
- AD: prefilter **32 DFC / 96 MTG**; literature candidates **183 DFC / 446 MTG** (union 529 unique).
- SCZ: prefilter **33**; literature candidates **118**.

**2. Deterministic PubMed pre-fetch (no LLM).** For every literature-track gene, an NCBI
E-utilities `esearch` retrieves the total disease hit count and the top-5 most-relevant PMIDs
for the disease query (2000–present), followed by a single batched `efetch` of those
abstracts. Because PMIDs and abstract text are pulled directly from NCBI, **hallucinated
citations are impossible by construction** — this replaces the separate PMID-verification
agent of the prior pipeline. Requests are paced at ~2.5/s with retry/back-off (no API key).
- AD: of 529 genes, 273 returned ≥1 abstract; 256 had zero literature and were auto-excluded with no model call.
- SCZ: of 118 genes, 42 returned ≥1 abstract; 76 had zero literature and were auto-excluded.

**3. LLM classification (Claude Sonnet 5, batched).** The pre-fetched abstracts are classified
in batches of ~20–30 genes. For each gene the model decides whether ≥1 provided abstract
establishes a **direct** link to the disease via one of its mechanism categories. Co-mention,
"expressed in brain", differential expression without a disease-specific finding,
other-disease context, and gene-symbol homonyms are excluded; supporting PMIDs are chosen only
from the provided real set.
- **AD — 13 categories:** APP processing/trafficking; Aβ generation/aggregation/clearance; tau
  phosphorylation/clearance; autophagy-lysosome; ubiquitin-proteasome system; mitochondrial
  dysfunction; endosomal-lysosomal trafficking; ER stress/UPR/ERAD; axonal transport; synaptic
  failure/loss; neuroinflammation; Ca²⁺ dysregulation; AD GWAS. Homonyms excluded include
  ATR→"ATR-FTIR", SACS→"single-atom catalysts", NNT→"number needed to treat". Symbol pass:
  273 → **109 included**.
- **SCZ — 14 categories:** dopaminergic signaling; glutamatergic/NMDA hypofunction; GABAergic
  interneuron; synaptic vesicle cycling/neurotransmitter release; dendritic spine
  morphology/density; neurodevelopmental processes; myelination/oligodendrocyte;
  mitochondrial dysfunction; neuroinflammation/complement/immune; epigenetic regulation; RNA
  processing/splicing; SCZ GWAS (common-variant); calcium-channel signaling; rare coding/CNV
  variant burden. (The last two extend the original 12-category SCZ set.) A known false-positive
  symbol, GK ("gatekeeping"/author initials), is correctly excluded. Symbol pass: 42 → **18 included**.

**4. Alias augmentation.** The strict `"SYMBOL"[TIAB]` query misses genes whose disease
literature uses an alternate symbol. For each not-yet-included gene the query is re-run as
`("SYMBOL"[TIAB] OR "ALIAS"[TIAB] …) AND "<disease>"[TIAB]` using NCBI Gene "Also known as"
aliases; only abstracts new relative to the symbol query are fetched, and only the genes that
gained ≥1 new abstract are re-classified. An explicit **alias-specificity guard** rejects any
alias shared with a different entity.
- AD: added **24 genes → 133 literature-included** (recovering e.g. MAPK8→"JNK", PAFAH1B1→"LIS1";
  correctly rejecting DCTN4→"P62"=SQSTM1 and SACM1L→"SAC1"=SYNJ1). Every one of the 433
  literature references is a real NCBI-fetched PMID.
- SCZ: 18 not-yet-included genes gained alias abstracts, but **all 18 were rejected by the
  specificity guard** (0 added) — every candidate alias was a homonym for a different entity or
  a generic term (e.g. SLC38A1→"SAT1" = a distinct gene; IARS1→"IRS" = immune-response system;
  SLK→"FYN"; PIK3CB→"PI3K"; RPS29→a SNP inside GABRB2). The literature track stays at **18**.

**5. Curator preservation.** Even alias augmentation cannot reach genes whose disease literature
uses a **protein-family name that NCBI does not list as a gene alias**.
- AD: auditing newly-excluded genes against the prior human-curated database and PubMed
  protein-name counts identified **three** such genes — **PPP3CB** ("calcineurin" = 292 AD papers
  vs 0 for the symbol), **SPTAN1** ("spectrin" = 55, "fodrin" = 18 vs 2), **EDEM3** (ERAD,
  described generically). These family names are shared across paralogues and cannot be safely
  auto-attributed, so they are preserved from prior manual curation with their original PMIDs and
  flagged `[curator]` / "Curator-preserved". (DCTN4 and SACM1L were **not** preserved — their only
  hits were wrong-gene aliases.)
- SCZ: the same audit found **no** genes requiring preservation. The established SCZ genes in the
  overlap (e.g. RB1CC1, TCF4, SCN2A, CACNA2D1, DLG2, CNTNAP2) are all captured by the
  comprehensive prefilter database, and the genes dropped relative to the prior three-agent
  SCZ database were fully reachable and correctly excluded on the merits (drug-response-only,
  null findings, or co-mention). Curator count = **0**.

**6. Annotation.** Every final-table gene carries its total disease PubMed hit count (`esearch`
Count). AD additionally carries the Open Targets Platform AD association score (disease
MONDO_0004975, v4 GraphQL); this score is not computed for SCZ. Prefilter-track genes carry the
curated-database source as their reference.

**7. Output (minimal-file).**
- AD: `ad_db_summary_table_{dfc,mtg}.csv` (columns Gene, Region, N_papers, AD_mechanism,
  References, References_verified, AD_involvement_verified, Agent3_resolution,
  Pubmed_total_hits, OpenTargets_AD_score) + `ad_db.R` (`ad_db_dfc`, `ad_db_mtg` lists of
  `n` = paper count / `ref` = top PMID, grouped by mechanism). Final: **DFC 74** (39 literature +
  3 curator + 32 prefilter), **MTG 215** (116 literature + 96 prefilter + 3 curator).
- SCZ: `scz_db_summary_table_dfc.csv` (same schema minus `OpenTargets_*`; `SCZ_mechanism` /
  `SCZ_involvement_verified`) + `scz_db.R` (`scz_db_dfc`). Final: **DFC 51** (18 literature +
  33 prefilter; 0 curator).

Both databases feed their figure panel and a **Table S15** sheet (AD DFC/MTG, SCZ DFC).

**Reproducibility & limitations.** The PubMed retrieval, alias lookup, hit counts and Open
Targets scores are fully deterministic; the mechanism classification is LLM-based (Claude
Sonnet 5) and therefore not bit-reproducible across runs, though every retained citation is a
verifiable, NCBI-fetched PMID. Symbol- and alias-based queries structurally miss genes known
only by a shared protein-family name; this is mitigated — but not fully eliminated — by curator
preservation benchmarked against prior manual curation (needed for AD; not required for SCZ).

---

## Condensed (≈210 words)

Reproducible dCoPA genes — genes in modules significant across both single-nucleus datasets of
a disease (AD: Gabitto + Liu, DFC and MTG; SCZ: brainSCOPE CMC + SZBDMulti-Seq, DFC) — were
annotated for disease association by a single pipeline run separately for Alzheimer's disease
and schizophrenia. Genes present in a curated disease database (OMIM, Open Targets, ClinVar,
and, for AD, AlzForum) were included directly; all remaining genes entered a literature
pipeline. For each, NCBI E-utilities retrieved the top-5 most-relevant PubMed abstracts for
`"gene"[TIAB] AND "<disease>"[TIAB]` (2000–present); fetching real abstracts up front makes
citation hallucination impossible. Claude Sonnet 5 then classified each gene as disease-linked
only if a retrieved abstract established a direct mechanistic link to one of the disease's
mechanism categories (13 for AD, 14 for SCZ), rejecting co-mentions and gene-symbol homonyms and
citing only retrieved PMIDs. Genes missed by the exact symbol were re-queried with NCBI gene
aliases under a homonym-specificity guard (recovering, for AD, MAPK8→JNK and PAFAH1B1→LIS1;
adding none for SCZ, where every alias candidate was a homonym for a different entity). Three
AD genes whose literature uses a shared protein-family name (PPP3CB/calcineurin,
SPTAN1/spectrin, EDEM3/ERAD) were preserved from manual curation; SCZ needed none. Final
databases: AD DFC 74 / MTG 215; SCZ DFC 51.

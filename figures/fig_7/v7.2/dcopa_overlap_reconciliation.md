# dCoPA overlap reconciliation & Table S15 construction record

**Date:** 2026-06-25
**Author of analysis:** Claude Code session
**Scope:** Why 5 MTG genes lacked module/subclass annotation in Table S15, the version-mismatch
investigation it triggered, the byte-for-byte regeneration of the dCoPA overlap files that
resolved it, and the resulting decision to drop those 5 rows.

---

## 1. Table S15 — what it is

`table_s15/v1/format_table_s15.R` builds `Kang_Table_S15_v1.xlsx`:

- Combines the two AD literature-review summary tables into one workbook, **one tab per region**
  (MTG, DFC), with a **Legend** tab first.
  - Source: `fig_7/v7.2/ad_db_summary_table_{mtg,dfc}.csv`
- Adds dCoPA context columns: **Module index (CTRL)**, **Module index (AD)**,
  **Associated subclass**, **Direction**.
  - Source for these: `fig_7/v7.1/panel_B_dcopa_genelist.csv` restricted to the reproducible
    (Gene, Celltype) pairs in `fig_7/v7.2/{mtg,dfc}_overlaps.csv`.
- Direction relabelled to pathological framing ("Lower/Higher in pathological samples").
- `N_papers` column removed; `References` PMIDs spaced as `a | b | c`.

### N_papers / References provenance (recorded for completeness)
`N_papers` and `References` come from the three-agent literature pipeline (`ad_db.R` + prompts).
`N_papers` (agent-tallied supporting papers) and the number of listed `References` PMIDs **agree
exactly only for genes reviewed under the v3 prompt** (MTG, 82 genes; perfect 1:1, "5+" capped at
5). For genes reviewed under the v2 prompt (all DFC + the original ~63 MTG genes), `N_papers`
systematically exceeds the listed PMIDs (e.g. `N_papers=5` → 2 PMIDs). This is why `N_papers` was
dropped from S15.

---

## 2. The trigger: 5 MTG genes with no annotation

5 MTG summary genes had blank Module index / subclass / direction:
**CAMTA1, CAPZB, CHSY1, CLCN6, CLK1**. They are all confirmed *literature-pipeline* genes
(verified PMIDs, Agent3 = "Confirmed") but are **absent from `mtg_overlaps.csv`**, so the
annotation join returned nothing.

---

## 3. Version-mismatch investigation

| File | Role | Date |
|---|---|---|
| `v7.2/mtg_overlaps.csv` (== `v7.1/mtg_overlaps.csv`, md5 `70377c7c…`) | reproducibility set used by S15 | 2026-03-25 run |
| `v7.1/panel_B_dcopa_genelist.csv` | module annotation source | 2026-06-17 (regenerated; **added the `Module` column** absent in the v7 March panel_B) |

So the reproducibility filter (March) and the module annotation (June) came from **different
runs**. Against the June `panel_B`, CLCN6 *looked* reproducible (present in both Gabitto and Liu,
same modules CTRL 337 / AD 149, sharing subclasses L2/3 IT + Vip) — conflicting with its absence
from the overlap file.

### Why panel_B cannot decide reproducibility
The authoritative criterion (from `fig7_v7.1.R`) is much stricter than "gene appears in both
datasets":

1. Per dataset, load the per-subclass euclidean-distance dCoPA tables and filter to
   `mod %in% these_mods` (module size > 3 and bulk-correlation significant), `sig_FDR == TRUE`,
   and `Consistency %in% c(0, 1)` (fully consistent).
2. Cross-dataset overlap (`find_output_overlap`): inner-join **Gabitto** and **Liu** on
   `(mod, Direction, Celltype, Consistency)`.
3. Map surviving modules → genes (CTRL via `mods`, AD via `mods_case`) → the two
   `panel_C_D_*_dcopa_genelist.csv` (CTRL-module and AD-module overlaps).
4. `find_comps` inner-joins the CTRL-module and AD-module genelists on `(Celltype, Gene)`.

**Reproducible = a gene present, for the same subclass, in BOTH a reproducible CTRL module AND a
reproducible AD module**, each itself a `sig_FDR` + `Consistency`-filtered Gabitto∩Liu
intersection. `panel_B` has neither `sig_FDR` nor `Consistency`, so a panel_B-only reconstruction
was ~6× too large (3402 vs 544 genes) — not faithful.

---

## 4. Regeneration (the decisive test)

`regen_overlaps.R` (this folder) copies the overlap-generation chain **verbatim** from
`fig7_v7.1.R` (plotting/GSEA omitted). All upstream inputs were present and stable:

| Input | Date |
|---|---|
| CTRL kme table (`SEAAD2024_AllADVsCon_DFC`) | 2025-02-13 |
| AD kme table (`rosmap_AD_…`) | 2025-12-24 |
| 8 euclidean-distance tables | 2026-03-22 |
| `bulk_cors_sigcount_bonf_1158.csv` | 2026-03-03 |

### Result — regenerated vs committed overlap files

| File | Regenerated md5 | Committed `v7.2` md5 | Match | Rows |
|---|---|---|---|---|
| `mtg_overlaps.csv` | `70377c7c…` | `70377c7c…` | **identical** | 848 |
| `dfc_overlaps.csv` | `918875d1…` | `918875d1…` | **identical** | 254 |

Genes added: **0**. Genes removed: **0**. **CLCN6 and all 5 genes are absent from the regenerated
set**, exactly as in the committed file.

### Adjudication
- The committed `mtg_overlaps.csv` is **authoritative and current** — fully reproducible from the
  upstream data, **not stale**.
- All 5 genes are genuinely **non-reproducible** by the authoritative criterion.
- Nuance: "non-reproducible" ≠ "not a dCoPA gene in both datasets". CLK1 isn't a dCoPA gene at all;
  CAPZB is Liu-only; CHSY1 is Gabitto-only; **CAMTA1 and CLCN6 *are* dCoPA genes in both datasets**
  but fail the stricter `sig_FDR` + `Consistency` + CTRL∩AD-per-subclass test. The June `panel_B`
  "reproducible" appearance for CLCN6 was an artifact of that looser view.

| Gene | dCoPA in both datasets (panel_B)? | In reproducible overlap (authoritative)? |
|---|---|---|
| CLK1 | no (absent from panel_B) | no |
| CAPZB | Liu only | no |
| CHSY1 | Gabitto only | no |
| CAMTA1 | yes (same modules, different subclasses) | no |
| CLCN6 | yes (same modules, shares subclasses) | no |

---

## 5. Actions taken (2026-06-25)

1. **Dropped the 5 non-reproducible rows from Table S15.** `format_table_s15.R` now keeps only
   summary genes present in the overlap set (`s <- s[s$Gene %in% ann$Gene, ]`) and logs which
   genes are dropped. Rebuilt `Kang_Table_S15_v1.xlsx`:
   - MTG: 222 → **217** rows; DFC: **77** rows (unchanged). No blank annotation cells remain.
2. **Merged the verified regeneration into v7.2.** `regen_overlaps.R` was moved here
   (`fig_7/v7.2/`), making v7.2 self-contained for overlap generation (previously it only *read*
   the overlap files, which had been carried over from v7.1). The committed `v7.2/{mtg,dfc}_overlaps.csv`
   were left untouched because the regeneration reproduces them byte-for-byte.
3. **Removed the scratch folder** `fig_7/overlap_regen_test/` after merging its relevant parts
   (the regen script → `v7.2/`; this findings record). Its other contents were duplicates
   (overlap CSVs identical to `v7.2`) or regenerable intermediates (`panel_C_D_*` genelists).

## 6. How to reproduce
```
cd Code_for_figures/fig_7/v7.2
Rscript regen_overlaps.R          # writes mtg/dfc_overlaps.csv + panel_C_D_* genelists here
cd ../../table_s15/v1
Rscript format_table_s15.R        # rebuilds Kang_Table_S15_v1.xlsx (drops non-reproducible genes)
```

# AD/SCZ module-overlap p-value: the four runs

Summary of the four statistical approaches used to test whether Alzheimer's
disease (AD) and schizophrenia (SCZ) converge on the same dCoPA-perturbed modules
more than expected by chance. Runs 1–3 are one-sided Fisher's exact tests in
`ad_scz_comparison.R`; run 4 is the permutation test in
`ad_scz_permutation_test.R`. The four runs form a chain: each one's drawback
motivates the next.

**Shared inputs.** AD dCoPA = 63 modules / 124 (module, cell type, direction)
triplets. SCZ dCoPA = 65 modules / 162 triplets. Module universe *N* = 1016.
"Strict shared" = a (module, cell type, direction) triplet present in BOTH
diseases = **18 triplets**, coming from **12 unique modules** (all direction =
down). Module-level co-occurrence (a module hit in both diseases, any ct/dir) =
32 modules.

| # | Run | Unit of analysis | "Shared" count | Population | Odds ratio | p-value |
|---|---|---|---|---|---|---|
| 1 | Fisher, module membership | module | 32 | 1016 | 28.5 | 3.03e-25 |
| 2 | Fisher, strict triplet | (mod, ct, dir) triplet | 18 | 48,768 | 57.2 | 1.39e-24 |
| 3 | Fisher, de-pseudoreplicated (module level) | module | 12 | 1016 | 4.0 | 3.42e-04 |
| 4 | Permutation — Test 1 (global convergence) | module / triplet | 18 trip / 12 mod | 1016 (empirical null) | — | < 1e-5 |
| 4 | Permutation — Test 2 (conditional matching) | module | 12 of 32 | empirical null | — | 0.26 (n.s.) |

---

## Run 1 — Fisher's exact test, module membership

**What it does.** Treats each module as one draw. A module is a "success" for a
disease if it has any dCoPA hit in that disease (cell type / direction ignored).
2×2 table of {in AD, out} × {in SCZ, out} over the 1016-module universe: (32, 31,
33, 920). One-sided (`alternative = "greater"`) test for enrichment of joint hits.

**Rationale.** The simplest, most defensible-in-assumptions framing of the
question "do AD and SCZ hit the same *modules*?" Modules are the natural
independent unit (one label each), so Fisher's exchangeability assumption is least
strained here.

**Drawback → Run 2.** It ignores cell type and direction entirely, so it cannot
speak to the "shared" definition the downstream expression / PC1 analyses actually
use — a (module, cell type, direction) triplet matched across diseases. To test
enrichment in the space where the biology is defined, the next run moves the unit
of analysis down to the triplet.

## Run 2 — Fisher's exact test, strict triplet

**What it does.** Uses the strict definition of "shared" — a (module, cell type,
direction) triplet in both diseases (18) — as the drawn unit. Population is every
triplet that *could* be drawn: 1016 modules × 24 cell types × 2 directions =
48,768. Table (18, 106, 144, 48,500), one-sided.

**Rationale.** Matches the "shared" definition actually used to select modules for
the downstream analyses (`mod_group_defs.R`), so it formally tests enrichment in
the same space the biology is defined in.

**Drawback → Run 3.** It violates Fisher's independence assumption in two ways,
both of which make the p-value anti-conservative (too small):
1. **Pseudoreplication** — the 18 triplets come from only 12 modules, because a
   perturbed module tends to hit several cell types at once. They are ~12
   independent events counted as 18, which understates the variance of the overlap.
2. **Fictitious direction multiplier** — direction is near-deterministic (all
   shared triplets are "down"), so the ×2 directions in the 48,768 population is
   largely imaginary and artificially shrinks the expected-by-chance overlap.

Both biases pull the same way, so the p (~1e-24) should be treated as an
**upper-bound heuristic**, not the headline. Run 3 removes both by returning to the
module as the unit while keeping the matched definition.

## Run 3 — Fisher's exact test, de-pseudoreplicated (module level)

**What it does.** Keeps the strict cell-type+direction-*matched* definition of
"shared" but collapses back to the module as the unit: shared = a module with a
ct+dir-matched triplet in BOTH diseases = 12. Margins are module-level dCoPA
membership (63, 65) over 1016. Table (12, 51, 53, 900), one-sided. A module hit in
both diseases but never in a matched (ct, dir) lands off-diagonal — counted as NOT
shared (deliberately conservative).

**Rationale.** Directly repairs the two biases in Run 2: one vote per module
removes the pseudoreplication, and dropping the cell-type × direction population
multipliers removes the fictitious inflation. Expected overlap by chance rises
from ~0.4 (Run 2) to ~4 modules; observed = 12. The ~20-order-of-magnitude gap
between Run 2 and Run 3 (1e-24 → 3e-4) *is* the pseudoreplication + direction
inflation made concrete. This is the defensible closed-form p-value.

**Drawback → Run 4.** Even this repaired Fisher test still assumes modules are
equiprobable, exchangeable hits — a closed-form null that cannot encode the real
data structure (some modules are far more "hittable" than others; footprints
cluster). Crucially, it also **cannot separate two distinct signals**: (a) AD and
SCZ landing on the same modules at all, versus (b) them landing on the same
modules in the same *cell type and direction*. A permutation test can build the
null from the real footprints and, with two complementary designs, tease these
apart.

## Run 4 — Permutation test (empirical null)

The permutation test builds the null from the **real data structure** — it
preserves each disease's observed module clustering and ct/dir footprint and
randomises only the cross-disease association, so the dependence is preserved
rather than assumed away. B = 100,000, seed = 1. Two complementary nulls:

**Test 1 — global convergence.** Relabel SCZ's 65 modules to a random size-65
subset of the 1016 universe, carrying each module's real ct/dir footprint; AD held
fixed. Recount matched triplets/modules. *Rationale:* the rigorous replacement for
the Fisher enrichment tests — it asks whether the observed matched overlap beats
what you get when SCZ hits random modules while keeping its real footprint.
Result: obs 18 triplets / 12 modules vs null means ~1.2 / ~1.0, 0/100,000 ≥
observed, **p < 1e-5**. Confirms strong convergence; the rigorous null mean (~1)
is *below* the Fisher expected (~4), so if anything the closed-form tests
understated it.

**Test 2 — conditional matching.** Restrict to the 32 modules hit in BOTH
diseases; hold AD footprints fixed and permute which SCZ footprint pairs with which
module. Count modules whose AD/SCZ footprints share ≥1 (ct, dir). *Rationale:*
isolates cell-type/direction *specificity* from the module co-occurrence signal
that Test 1 and all the Fisher runs bundle together. Result: obs 12 of 32 vs null
mean ~10.3, **p = 0.26 (not significant)** — given that a module is hit in both
diseases, the ct/dir matching is no better than random pairing of the observed
footprints (hits are overwhelmingly "down" in a few overlapping cell types).

**Drawback (endpoint).** The permutation approach costs interpretability and
compute for its rigor: results are Monte-Carlo estimates (the reported p < 1e-5 is
really 1/(B+1) = 9.9999e-06, floored by B), and the two tests answer *different*
questions, so they must be cited carefully rather than as one number. Test 2's
null result is the terminal caveat of the whole chain: the ct/dir-resolved sharing
is real enough to *select* modules on, but is not independently significant.

---

## Bottom line (how to cite)

- **Report the module-level convergence** (Run 3 Fisher, p ≈ 3e-4; corroborated by
  permutation Run 4 / Test 1, p < 1e-5) as the enrichment result.
- The strict-triplet Fisher p (~1e-24, Run 2) is inflated by pseudoreplication + an
  oversized population — use only as an upper-bound heuristic.
- The cell-type/direction-resolved "shared" definition is a valid way to **select**
  the module set for the downstream expression / PC1 analyses, but permutation
  Run 4 / Test 2 shows it is **not** independently significant evidence of
  ct/dir-specific sharing — describe that pattern descriptively, not as a finding.

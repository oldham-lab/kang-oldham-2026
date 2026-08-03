# Permutation tests for AD/SCZ module overlap under the STRICT "shared"
# definition (a (module, cell type, direction) triplet present in both diseases;
# see mod_group_defs.R and ad_scz_comparison.R).
#
# WHY A PERMUTATION TEST: the closed-form Fisher tests in ad_scz_comparison.R
# assume the drawn items are independent, exchangeable, equiprobable units. The
# triplet test violates this (triplets cluster by module; direction is near-
# deterministic), and even the module-level test still assumes modules are
# equiprobable hits and cannot isolate cell-type/direction matching from mere
# module co-occurrence. A permutation test sidesteps all of this by building the
# null from the REAL data structure: we keep each disease's observed footprint
# (its module-level clustering and ct/dir pattern) intact and only randomise the
# cross-disease association, so the dependence structure is preserved rather than
# assumed away.
#
# Two complementary nulls:
#   TEST 1 (global convergence): does AD/SCZ matched overlap exceed what is seen
#     when SCZ hits a RANDOM set of modules (carrying its real ct/dir footprints)
#     independent of AD? Rigorous replacement for the Fisher enrichment tests.
#   TEST 2 (conditional matching): among modules hit in BOTH diseases, do AD and
#     SCZ align on the SAME (ct, dir) more than expected if their footprints were
#     randomly paired across those modules? Isolates ct/dir specificity from the
#     module co-occurrence signal that TEST 1 (and the Fisher tests) bundle in.

library(tidyverse)
library(data.table)

set.seed(1)
B <- 100000   # permutations per test

# --- Data ------------------------------------------------------------------
ad  <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/ad_dfc.csv"),  data.table = F)
scz <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/scz_dfc.csv"), data.table = F)

# 1016-module universe (identical construction to mod_group_defs.R).
mdir   <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme <- fread(file.path(mdir, "kme_tables", "topmodposbc_table.csv"), data.table = F)
mods   <- tapply(datkme[[2]], datkme[[3]], list)
filter_under <- 3
these_mods <- as.numeric(names(mods)[sapply(mods, length) > filter_under])
sigcount   <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"), data.table = F)
these_mods <- these_mods[!these_mods %in% which(sigcount$vals < 2)]
stopifnot(length(these_mods) == 1016)

# Distinct (mod, ct, dir) triplets per disease, restricted to the universe.
adk  <- ad  %>% filter(mod %in% these_mods) %>% distinct(mod, Celltype, Direction)
sczk <- scz %>% filter(mod %in% these_mods) %>% distinct(mod, Celltype, Direction)

# Per-module footprint = the set of "ct|dir" tags a module has a dCoPA hit in.
tag      <- function(d) paste(d$Celltype, d$Direction, sep = "|")
ad_fp    <- split(tag(adk),  adk$mod)    # named by module id (character)
scz_fp   <- split(tag(sczk), sczk$mod)
ad_ids   <- names(ad_fp)                 # 63
scz_ids  <- names(scz_fp)                # 65

# --- Observed statistics ---------------------------------------------------
matched_per_scz <- function(ids, fps) {  # ids: assigned module id per scz module
  trip <- 0L; modu <- 0L
  for (i in seq_along(fps)) {
    af <- ad_fp[[ ids[i] ]]
    if (!is.null(af)) {
      k <- sum(fps[[i]] %in% af)
      if (k > 0) { trip <- trip + k; modu <- modu + 1L }
    }
  }
  c(trip = trip, modu = modu)
}
obs1 <- matched_per_scz(scz_ids, scz_fp)
obs_trip <- obs1["trip"]   # 18 shared triplets
obs_mod  <- obs1["modu"]   # 12 shared modules
cat(sprintf("Observed: %d matched triplets, %d matched modules\n\n", obs_trip, obs_mod))

# ===========================================================================
# TEST 1 -- global convergence null
# Relabel SCZ's 65 modules to a random size-65 subset of the 1016 universe,
# carrying each module's real ct/dir footprint. AD held fixed. This preserves
# SCZ's within-module clustering and direction structure exactly; it breaks only
# the cross-disease module identity. Recount matched triplets and modules.
# ===========================================================================
null_trip <- integer(B); null_mod <- integer(B)
universe  <- as.character(these_mods)
for (b in 1:B) {
  newids <- sample(universe, length(scz_fp))      # 65 distinct ids
  m <- matched_per_scz(newids, scz_fp)
  null_trip[b] <- m["trip"]; null_mod[b] <- m["modu"]
}
p_trip <- (sum(null_trip >= obs_trip) + 1) / (B + 1)
p_mod  <- (sum(null_mod  >= obs_mod ) + 1) / (B + 1)

cat("TEST 1 -- global convergence (SCZ modules randomly relocated)\n")
cat(sprintf("  matched triplets: obs=%d  null mean=%.3f  null max=%d  n>=obs=%d  p=%.6g\n",
            obs_trip, mean(null_trip), max(null_trip), sum(null_trip >= obs_trip), p_trip))
cat(sprintf("  matched modules : obs=%d  null mean=%.3f  null max=%d  n>=obs=%d  p=%.6g\n\n",
            obs_mod, mean(null_mod), max(null_mod), sum(null_mod >= obs_mod), p_mod))

# ===========================================================================
# TEST 2 -- conditional matching null
# Restrict to modules with a dCoPA hit in BOTH diseases. Hold AD footprints fixed; permute
# which SCZ footprint is paired with which module. Count modules whose AD and
# SCZ footprints share >=1 (ct, dir). If observed >> null, the matching is
# module-specific (same module hit the same way); if observed ~ null, the matches
# are just what random pairing of these footprints would give (e.g. because most
# hits are "down" in a few cell types). This isolates ct/dir specificity from the
# module co-occurrence enrichment.
# ===========================================================================
co_ids   <- intersect(ad_ids, scz_ids)            # 32 modules with a dCoPA hit in both
ad_co    <- ad_fp[co_ids]
scz_co   <- scz_fp[co_ids]
n_co     <- length(co_ids)
overlap  <- function(a, s) any(s %in% a)
obs2     <- sum(mapply(overlap, ad_co, scz_co))   # should equal obs_mod (12)

null2 <- integer(B)
for (b in 1:B) {
  perm <- sample(n_co)
  null2[b] <- sum(mapply(overlap, ad_co, scz_co[perm]))
}
p2 <- (sum(null2 >= obs2) + 1) / (B + 1)

cat("TEST 2 -- conditional matching (SCZ footprints re-paired among the 32 co-hit modules)\n")
cat(sprintf("  matched modules : obs=%d of %d co-hit  null mean=%.3f  null max=%d  n>=obs=%d  p=%.6g\n",
            obs2, n_co, mean(null2), max(null2), sum(null2 >= obs2), p2))

# ===========================================================================
# RESULTS (seed = 1, B = 100,000)
# Observed: 18 matched triplets, 12 matched modules.
#
# TEST 1 (global convergence):
#   matched triplets: obs=18  null mean=1.235  max=11  p < 1e-5
#   matched modules : obs=12  null mean=0.996  max=8   p < 1e-5
#   => AD/SCZ converge on the same modules FAR beyond chance, even with each
#      disease's module clustering and ct/dir footprints held intact. This is the
#      defensible confirmation of the enrichment; note the rigorous null mean (~1)
#      is below the Fisher table's expected (~4), so the closed-form tests if
#      anything understated the convergence.
#
# TEST 2 (conditional matching):
#   matched modules : obs=12 of 32 co-hit  null mean=10.310  max=18  p=0.26
#   => NOT significant. Given that a module is hit in both diseases, the ct/dir
#      matching (12/32) is no more than random pairing of the observed footprints
#      yields (~10.3). The matching is explained by marginal structure (hits are
#      overwhelmingly "down" in a few overlapping cell types), NOT by module-
#      specific concordance.
#
# INTERPRETATION: cite the MODULE-LEVEL convergence (Test 1) as the result. The
# cell-type/direction-resolved "shared" definition is a valid way to SELECT the
# module set for the downstream expression / PC1 analyses, but it is NOT evidence
# of ct/dir-specific sharing -- describe that pattern descriptively, not as an
# independently significant finding. The strict-triplet Fisher p (~1e-24) was
# inflated by pseudoreplication + an oversized population, not real specificity.
# ===========================================================================

# ---------------------------------------------------------------------------
# TECHNICAL SUMMARY OF THE FOUR RUNS
# ---------------------------------------------------------------------------
# Inputs: AD dCoPA = 63 modules / 124 (module, cell type, direction) triplets;
# SCZ dCoPA = 65 modules / 162 triplets; module universe N = 1016. Strict "shared"
# (triplet matched on module, cell type, and direction) = 18 triplets / 12 modules.
# Runs 1-3 are one-sided (greater) Fisher's exact tests (ad_scz_comparison.R);
# run 4 is permutation-based (this script, B = 100,000, seed = 1).
#
# Run 1 -- Fisher, module membership. Unit = module; shared = module with a dCoPA hit in both
#   diseases (any cell type/direction) = 32. Table (32, 31, 33, 920), N = 1016.
#   Odds ratio = 28.47.  p = 3.031297e-25.
#
# Run 2 -- Fisher, strict triplet. Unit = (module, cell type, direction) triplet;
#   shared = 18. Population = 1016 modules x 24 cell types x 2 directions = 48,768.
#   Table (18, 106, 144, 48500).  Odds ratio = 57.16.  p = 1.391763e-24.
#
# Run 3 -- Fisher, de-pseudoreplicated (module level). Unit = module; shared =
#   module with a cell-type+direction-matched triplet in both = 12; margins =
#   module-level dCoPA membership (63, 65), N = 1016. Table (12, 51, 53, 900).
#   Odds ratio = 3.987.  p = 3.422739e-04.
#
# Run 4 -- Permutation tests (B = 100,000, seed = 1).
#   Test 1, global convergence: SCZ's 65 modules relabeled to a random size-65
#     subset of the 1016 universe with footprints attached, AD fixed.
#     Matched triplets: obs = 18, null mean = 1.235, null max = 11, 0/100,000
#       >= obs, p = 9.9999e-06.
#     Matched modules:  obs = 12, null mean = 0.996, null max = 8, 0/100,000
#       >= obs, p = 9.9999e-06.
#   Test 2, conditional matching: among the 32 modules with a dCoPA hit in both, SCZ footprints
#     permuted across modules, AD fixed. Matched modules: obs = 12, null mean =
#     10.310, null max = 18, 25,881/100,000 >= obs, p = 0.258817.
# ---------------------------------------------------------------------------

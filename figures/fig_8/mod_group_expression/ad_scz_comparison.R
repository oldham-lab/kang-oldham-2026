library(tidyverse)
library(data.table)

ad <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/ad_dfc.csv"))
scz <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/scz_dfc.csv"))

length(unique(ad$mod))
# 63
length(unique(scz$mod))
# 65

length(intersect(unique(ad$mod), unique(scz$mod)))
# 32

# Mods shared by AD and SCZ on the (mod, Celltype, Direction) key.
ad_keys  <- ad  %>% distinct(mod, Celltype, Direction)
scz_keys <- scz %>% distinct(mod, Celltype, Direction)

intersect_tbl <- inner_join(ad_keys, scz_keys,
                            by = c("mod", "Celltype", "Direction")) %>%
  arrange(mod, Celltype, Direction)

intersect_tbl
nrow(intersect_tbl)
# 18

write_csv(intersect_tbl,
          file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/ad_scz_mod_celltype_direction_intersect.csv"))

# --- Fisher's exact test for overlap of AD and SCZ module sets ---
# Population: 1016 modules (labeled 1..1016) from which each disease draws.
N_total <- 1016

ad_mods  <- unique(ad$mod)
scz_mods <- unique(scz$mod)

n_ad   <- length(ad_mods)                       # 63
n_scz  <- length(scz_mods)                       # 65
n_both <- length(intersect(ad_mods, scz_mods))   # 32

# 2x2 contingency table:
#                in SCZ      not in SCZ
#   in AD        n_both      n_ad  - n_both
#   not in AD    n_scz-both  N - (union)
ctab <- matrix(c(n_both,
                 n_ad  - n_both,
                 n_scz - n_both,
                 N_total - (n_ad + n_scz - n_both)),
               nrow = 2, byrow = TRUE,
               dimnames = list(AD  = c("in", "out"),
                               SCZ = c("in", "out")))
ctab

fisher.test(ctab, alternative = "greater")$p.value

# Fisher's exact test (one-sided, testing for enrichment): p < 2.2e-16, odds ratio â 28.5, 95% CI [16.5, Inf).
# exact p-value: 3.031297e-25

# --- Fisher's exact test, STRICT "shared" definition ---------------------
# Same enrichment test, but an item is now a (module, cell type, direction)
# TRIPLET rather than a bare module -- i.e. the definition of "shared" used in
# mod_group_expression/ and pc1_variance/ (see mod_group_defs.R): a triplet
# present in BOTH AD and SCZ (18 such triplets, all direction = down).
#
# Population: every triplet that could possibly be drawn =
#   1016 modules x 24 cell types x 2 directions = 48,768.
# (24 = the full cell-type panel over which the dCoPA test was run; every tested
#  cell type is a slot a shared triplet could have landed in, including ones
#  with zero hits. 2 = directions +1 / -1. The result is essentially invariant
#  to this count: 21/23/24 cell types all give p ~ 1e-23 to 1e-24.)
N_celltypes  <- 24
N_dir        <- 2
N_total_trip <- N_total * N_celltypes * N_dir      # 48,768

n_ad_trip   <- nrow(ad_keys)                        # 124 distinct AD triplets
n_scz_trip  <- nrow(scz_keys)                       # 162 distinct SCZ triplets
n_both_trip <- nrow(intersect_tbl)                  # 18 shared triplets

ctab_trip <- matrix(c(n_both_trip,
                      n_ad_trip  - n_both_trip,
                      n_scz_trip - n_both_trip,
                      N_total_trip - (n_ad_trip + n_scz_trip - n_both_trip)),
                    nrow = 2, byrow = TRUE,
                    dimnames = list(AD  = c("in", "out"),
                                    SCZ = c("in", "out")))
ctab_trip

fisher.test(ctab_trip, alternative = "greater")$p.value

# Strict Fisher's exact test (one-sided, enrichment): odds ratio ~ 57.2.
# expected overlap by chance ~ 0.41 triplets; observed = 18.
# exact p-value: 1.391763e-24

# --- Fisher's exact test, DE-PSEUDOREPLICATED (module-level) --------------
# Rationale: the strict (triplet) test above violates Fisher's core assumption
# that the items being drawn are INDEPENDENT, exchangeable units. They are not:
#   (1) A perturbed module tends to have a dCoPA hit across several cell types at once, so
#       its triplets co-occur. The 18 shared triplets come from only 12 unique
#       modules -- they are ~12 independent events, not 18. Counting them as 18
#       independent successes understates the variance of the overlap.
#   (2) Direction is near-deterministic given (mod, ct) (all shared triplets are
#       down), so the x2 direction factor in the triplet population is largely
#       fictitious and shrinks the expected-by-chance overlap artificially.
# Both biases push the triplet p-value the SAME way -- too small (anti-
# conservative). This run removes both by making the MODULE the unit of analysis
# (each module counts once) while keeping the strict cell-type+direction-matched
# definition of "shared":
#   in AD     = modules with a dCoPA hit in AD (any ct/dir)   -> 63
#   in SCZ    = modules with a dCoPA hit in SCZ (any ct/dir)  -> 65
#   shared    = modules with a ct+dir-MATCHED triplet in BOTH diseases -> 12
#   population = 1016 modules (no cell-type / direction multipliers)
# A module with a dCoPA hit in both diseases but never in a matched (ct, dir) lands in
# an off-diagonal cell, i.e. it is counted as NOT shared -- this is deliberately
# conservative: only matched modules score as successes. Expected overlap under
# the module-level margins is ~4, so 12 is still an enrichment, but far less
# extreme than the triplet test (p ~ 3e-4 vs ~1e-24). The ~20-order gap between
# the two is exactly the pseudoreplication + direction inflation, made concrete;
# treat this module-level p-value as the defensible one and the triplet p-value
# as an upper-bound heuristic. (For a fully rigorous null that keeps the within-
# disease clustering intact, permute module labels -- see prior discussion.)
N_total_mod <- N_total                              # 1016 modules

n_ad_mod   <- length(ad_mods)                       # 63 modules with a dCoPA hit in AD
n_scz_mod  <- length(scz_mods)                      # 65 modules with a dCoPA hit in SCZ
n_both_mod <- n_distinct(intersect_tbl$mod)         # 12 strict-shared modules

ctab_mod <- matrix(c(n_both_mod,
                     n_ad_mod  - n_both_mod,
                     n_scz_mod - n_both_mod,
                     N_total_mod - (n_ad_mod + n_scz_mod - n_both_mod)),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(AD  = c("in", "out"),
                                   SCZ = c("in", "out")))
ctab_mod

fisher.test(ctab_mod, alternative = "greater")$p.value

# De-pseudoreplicated Fisher's exact test (one-sided, enrichment): odds ratio ~ 4.0.
# expected overlap by chance ~ 4.03 modules; observed = 12.
# exact p-value: 0.0003422739


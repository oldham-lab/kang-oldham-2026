# Shared definitions for the fig_8 module-group analyses (mean_expression/ and
# pc1_variance/). Provides the Shared/Unshared/Neither grouping, the module
# universe + gene map, plotting helpers, and the raw bulk expression matrix.
# Analysis-specific data transforms live in each subfolder's own *_defs.R.
#
# "Shared" is defined per (module, cell type, direction): a (mod, ct, dir)
# triplet present in BOTH the AD and SCZ DE tables (18 such triplets, all
# direction = -1). Grouping is therefore CELL-TYPE-RESOLVED:
#   within a given cell type c,
#     Shared    = module is a shared triplet in c
#     Unshared  = module is DE in c in only one disease (or opposite directions)
#     Neither   = module (in the 1016 universe) is not DE in c
# A module can be Shared in one cell type and Neither in another.
#
# Sourcing this file provides:
#   gene_tbl            gene -> module (restricted to the 1016-module universe)
#   these_mods          the 1016-module universe
#   shared_keys         the 18 shared (mod, ct, Direction) triplets
#   shared_celltypes    cell types that contain >=1 shared triplet
#   expr_celltypes      the 23 cell-type labels
#   group_in_celltype(mod, ct)  -> "Shared"/"Unshared"/"Neither"/NA
#   bulk_expr           raw bulk expression matrix (gene x sample; ComBat-log)
#   sig_stars / signif_brackets / delta_heatmap   plotting helpers

library(tidyverse)
library(data.table)
library(ggsignif)

REGION <- "PFC"

# Significance stars from a (BH-adjusted) p-value.
sig_stars <- function(p) ifelse(p < 0.001, "***",
                         ifelse(p < 0.01,  "**",
                         ifelse(p < 0.05,  "*", "ns")))

# Build a bracket-annotation table for geom_signif from Wilcoxon results.
# `glevels` is the ordered x-axis group factor; xmin/xmax are the box positions.
# All other columns (e.g. the facet variable) are passed through unchanged.
signif_brackets <- function(wres, glevels) {
  wres %>%
    mutate(xmin = pmin(match(group1, glevels), match(group2, glevels)),
           xmax = pmax(match(group1, glevels), match(group2, glevels)),
           span = xmax - xmin,
           label = sprintf("%s\np=%s", sig_stars(p_adj_BH),
                           formatC(p_adj_BH, format = "g", digits = 2)))
}

# Lay out geom_signif brackets ABOVE each facet's boxplot whiskers, stacked with
# room for the (2-line) labels, so nothing overlaps the boxes/whiskers or each
# other. `d` = plotting data (needs `group` + the value column); `sig` from
# signif_brackets(). Returns list(sig = sig + y_position, ylim = c(lo, hi)).
bracket_layout <- function(sig, d, facet_col, value_col, group_col = "group", frac = 0.24) {
  whisk <- function(v, upper) {                       # boxplot whisker (1.5*IQR)
    q <- quantile(v, c(.25, .75), names = FALSE, na.rm = TRUE); iqr <- q[2] - q[1]
    if (upper) max(v[v <= q[2] + 1.5 * iqr], na.rm = TRUE)
    else       min(v[v >= q[1] - 1.5 * iqr], na.rm = TRUE)
  }
  ext <- d %>%
    group_by(.data[[facet_col]], .data[[group_col]]) %>%
    summarise(t = whisk(.data[[value_col]], TRUE),
              b = whisk(.data[[value_col]], FALSE), .groups = "drop_last") %>%
    summarise(top = max(t), bot = min(b), .groups = "drop")
  step <- (max(ext$top) - min(ext$bot)) * frac    # gap between stacked brackets (fits a 2-line label)
  sig <- sig %>%
    left_join(ext, by = facet_col) %>%
    group_by(.data[[facet_col]]) %>%
    arrange(span, .by_group = TRUE) %>%
    mutate(y_position = top + step * (row_number() - 0.6)) %>%   # 0.4*step margin above whisker
    ungroup()
  list(sig = sig, ylim = c(min(ext$bot), max(sig$y_position) + step * 1.3))
}

# Companion effect-size heatmap for the single-nucleus figures (where 4 dodged
# series per group make in-plot brackets ambiguous). Rows = `yvar`, x = series,
# facet = pairwise comparison; tile = Cliff's delta, text = BH-adj p stars.
delta_heatmap <- function(wres, yvar, ylevels = NULL, comp_levels = NULL,
                          title = NULL, caption = NULL) {
  d <- wres %>%
    mutate(series     = paste(dataset, condition),
           comparison = paste(group1, "vs", group2),
           star       = sig_stars(p_adj_BH))
  if (!is.null(comp_levels)) d$comparison <- factor(d$comparison, levels = comp_levels)
  if (!is.null(ylevels))     d[[yvar]]    <- factor(d[[yvar]], levels = ylevels)
  ggplot(d, aes(x = series, y = .data[[yvar]], fill = cliffs_delta)) +
    geom_tile(color = "grey85") +
    geom_text(aes(label = star), size = 4) +
    facet_wrap(~ comparison) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-1, 1), name = "Cliff's delta") +
    labs(x = NULL, y = NULL, title = title, caption = caption) +
    theme_bw(base_size = 16) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y  = element_text(size = 14),
          strip.text   = element_text(size = 15),
          plot.title   = element_text(size = 18, face = "bold"),
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 13),
          legend.position = "right")
}

# ---------------------------------------------------------------------------
# Module universe (1016) and gene -> module map
# ---------------------------------------------------------------------------
ad  <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/ad_dfc.csv"),  data.table = F)
scz <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/scz_dfc.csv"), data.table = F)

mdir     <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme   <- fread(file.path(mdir, "kme_tables", "topmodposbc_table.csv"), data.table = F)
gene_sym <- datkme[[2]]
gene_mod <- datkme[[3]]

mods         <- tapply(gene_sym, gene_mod, list)
modulelength <- sapply(mods, length)
filter_under <- 3
these_mods   <- as.numeric(names(mods)[modulelength > filter_under])
sigcount     <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"), data.table = F)
these_mods   <- these_mods[!these_mods %in% which(sigcount$vals < 2)]
stopifnot(length(these_mods) == 1016)

gene_tbl <- tibble(gene = gene_sym, mod = gene_mod) %>%
  filter(mod %in% these_mods) %>%
  distinct(gene, .keep_all = TRUE)

# ---------------------------------------------------------------------------
# Cell-type-resolved Shared / Unshared / Neither
# ---------------------------------------------------------------------------
# AD/SCZ celltype names differ from the expression columns only in case
# (e.g. PVALB -> Pvalb); map them onto the expression vocabulary.
expr_celltypes <- c("Astrocyte","Chandelier","Endothelial","L2/3 IT","L4 IT",
                    "L5 ET","L5 IT","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b",
                    "Lamp5","Lamp5 Lhx6","Microglia-PVM","OPC","Oligodendrocyte",
                    "Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip")
map_ct <- function(x) expr_celltypes[match(tolower(x), tolower(expr_celltypes))]

ad$ct  <- map_ct(ad$Celltype)
scz$ct <- map_ct(scz$Celltype)
stopifnot(!any(is.na(ad$ct)), !any(is.na(scz$ct)))

shared_keys <- inner_join(distinct(ad,  mod, ct, Direction),
                          distinct(scz, mod, ct, Direction),
                          by = c("mod", "ct", "Direction"))
stopifnot(nrow(shared_keys) == 18)
shared_celltypes <- sort(unique(shared_keys$ct))

shared_modct <- paste(shared_keys$mod, shared_keys$ct)
either_modct <- union(paste(ad$mod, ad$ct), paste(scz$mod, scz$ct))

# Vectorised group label for (mod, ct). NA if module is outside the 1016 universe.
group_in_celltype <- function(mod, ct) {
  key <- paste(mod, ct)
  out <- ifelse(key %in% shared_modct, "Shared",
         ifelse(key %in% either_modct, "Unshared", "Neither"))
  out[!mod %in% these_mods] <- NA
  factor(out, levels = c("Shared", "Unshared", "Neither"))
}

# Cell-type-AGNOSTIC ("overall") group, for collapsing the cell-type dimension:
# a module is Shared if it is a shared (mod, ct, dir) triplet in ANY cell type,
# Unshared if DE in ANY cell type in only one disease (or opposite directions),
# else Neither. Used by the *_allct figures that pool over all cell types.
# 12 Shared / 84 Unshared / 920 Neither modules. NA if outside the 1016 universe.
shared_mods_all <- sort(unique(shared_keys$mod))
either_mods_all <- union(ad$mod, scz$mod)
group_overall <- function(mod) {
  out <- ifelse(mod %in% shared_mods_all, "Shared",
         ifelse(mod %in% either_mods_all, "Unshared", "Neither"))
  out[!mod %in% these_mods] <- NA
  factor(out, levels = c("Shared", "Unshared", "Neither"))
}

# ---------------------------------------------------------------------------
# Raw bulk expression matrix (gene x sample; ComBat already log-normalized).
# Shared by both analyses: mean_expression/ takes per-gene row means;
# pc1_variance/ runs per-module PCA over the samples. Cols 1-2 are id + gene.
# ---------------------------------------------------------------------------
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table = F)

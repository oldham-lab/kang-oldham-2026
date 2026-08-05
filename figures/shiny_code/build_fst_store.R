# build_fst_store.R
# ---------------------------------------------------------------------------
# Final packaging step for the CoPA Shiny app's on-disk data store.
#
# The app (app.R) NO LONGER loads the big gene x cell-type matrices or the bulk
# correlation list into RAM. Instead they live on disk as `fst` files under
# `www/fst/` and are queried per-request (see `.fst_rows()` / `.bulk_module()`
# in app.R). This slashes startup RAM (~650 MB -> ~130 MB) and startup time
# (~6.6 s -> ~1.5 s), and removed the shinyapps.io out-of-memory crashes.
#
# Run this AFTER the upstream steps have assembled the source files under www/:
#   - make_shiny_snapshots_v1.3.R  -> www/CoPA_files/bulk_cor/bulk_cor_list.qs
#   - gene-projection pipeline     -> www/gene_projection/<ds>/PFC/{means,var,pct,n}/*.csv
#   - core_gbmap pipeline          -> www/core_gbmap/annotation_level_<1-4>/{means,var,n}/*.csv
#
# It produces (the ONLY data the app actually loads for these):
#   www/fst/gp/<ds>_{means,var,pct}.fst        (gene projection, per dataset)
#   www/fst/gb/annotation_level_<n>_{means,var}.fst   (core GBmap, per level)
#   www/fst/bulk_cor.fst                        (all modules stacked, row-indexed)
#   www/fst/store_meta.qs                       (tiny metadata: gene lists, cell
#     types, per-subclass n and global means, precomputed native quantiles,
#     bulk row index, and RELATIVE fst paths)
#
# Once built, the source CSVs / bulk_cor_list.qs are redundant for the deployed
# app and can be dropped from the bundle (keep them in staging for rebuilds).
# ---------------------------------------------------------------------------

suppressMessages({library(fst); library(qs); library(data.table)})

APP     <- Sys.getenv("SHINYAPP_DIR", "/home/gugene/ShinyApps/copacabana")           # app root (app's working dir)
GP_DS   <- c("SEAAD2024", "MIT", "Morabito2021", "LeinA9")
GB_LVL  <- paste0("annotation_level_", 1:4)

# metadata stores RELATIVE paths (app cwd == APP); we write to absolute paths.
rel <- function(...) file.path("www", "fst", ...)
abs <- function(relpath) file.path(APP, relpath)
dir.create(abs(rel("gp")), recursive = TRUE, showWarnings = FALSE)
dir.create(abs(rel("gb")), recursive = TRUE, showWarnings = FALSE)

# read a genes x cell-types CSV (first column = gene) into a matrix
read_mat <- function(p) {
  # header=TRUE is REQUIRED: some annotation levels have all-integer cell-type
  # labels (e.g. prelimCluster = 1..59), and without it fread mis-detects the
  # numeric header row as data, shifting genes and clobbering column names.
  d <- fread(p, data.table = FALSE, check.names = FALSE, header = TRUE)
  rn <- as.character(d[[1]]); d[[1]] <- NULL; rownames(d) <- rn
  as.matrix(d)
}
# write a matrix to fst with a leading `.gene` column (queried by .fst_rows)
write_mat_fst <- function(m, relpath) {
  df <- data.frame(.gene = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
  write_fst(df, abs(relpath), compress = 80)
}

meta <- list(gp = list(), gb = list())

## ---- gene projection (means / var / pct) ----------------------------------
gp_all_means <- c()
for (ds in GP_DS) {
  base <- file.path(APP, "www", "gene_projection", ds, "PFC")
  mn <- read_mat(file.path(base, "means", "genomewide_means.csv"))
  vr <- read_mat(file.path(base, "var",   "genomewide_var.csv"))
  pc <- tryCatch(read_mat(file.path(base, "pct", "genomewide_pct.csv")), error = function(e) NULL)
  nd <- fread(file.path(base, "n", "celltype_n.csv"), data.table = FALSE, check.names = FALSE)

  pm <- rel("gp", paste0(ds, "_means.fst")); write_mat_fst(mn, pm)
  pv <- rel("gp", paste0(ds, "_var.fst"));   write_mat_fst(vr, pv)
  pp <- NA_character_
  if (!is.null(pc)) { pp <- rel("gp", paste0(ds, "_pct.fst")); write_mat_fst(pc, pp) }

  meta$gp[[ds]] <- list(genes = rownames(mn), celltypes = colnames(mn),
                        n_vec = setNames(nd$n, nd$celltype),
                        global_mean = colMeans(mn, na.rm = TRUE),
                        means_path = pm, vars_path = pv, pct_path = pp)
  gp_all_means <- c(gp_all_means, as.numeric(mn))
  cat("gp", ds, "->", nrow(mn), "genes x", ncol(mn), "subclasses\n")
}
# union gene list + precomputed native-value quantiles (all genes) for y-axis coloring
meta$gp_genes <- sort(unique(unlist(lapply(meta$gp, function(x) x$genes))))
gv <- gp_all_means[is.finite(gp_all_means)]
meta$gp_native_qtl <- as.numeric(quantile(gv, probs = seq(0, 1, length.out = 101), na.rm = TRUE))
rm(gp_all_means, gv); gc()

## ---- core GBmap (means / var) ---------------------------------------------
for (lv in GB_LVL) {
  base <- file.path(APP, "www", "core_gbmap", lv)
  mn <- read_mat(file.path(base, "means", "genomewide_means.csv"))
  vr <- read_mat(file.path(base, "var",   "genomewide_var.csv"))
  nd <- fread(file.path(base, "n", "celltype_n.csv"), data.table = FALSE, check.names = FALSE)
  pm <- rel("gb", paste0(lv, "_means.fst")); write_mat_fst(mn, pm)
  pv <- rel("gb", paste0(lv, "_var.fst"));   write_mat_fst(vr, pv)
  meta$gb[[lv]] <- list(genes = rownames(mn), celltypes = colnames(mn),
                        n_vec = setNames(nd$n, nd$celltype),
                        global_mean = colMeans(mn, na.rm = TRUE),
                        means_path = pm, vars_path = pv)
  cat("gb", lv, "->", nrow(mn), "genes\n")
}

## ---- bulk_cor_list -> one stacked fst + per-module row index --------------
b <- qread(file.path(APP, "www", "CoPA_files", "bulk_cor", "bulk_cor_list.qs"), nthreads = 4)
nr    <- vapply(b, nrow, 0L)
ends  <- cumsum(nr); starts <- c(1L, head(ends, -1) + 1L)
write_fst(as.data.frame(data.table::rbindlist(b)), abs(rel("bulk_cor.fst")), compress = 80)
meta$bulk_path  <- rel("bulk_cor.fst")
meta$bulk_start <- as.integer(starts)
meta$bulk_end   <- as.integer(ends)
rm(b); gc()
cat("bulk_cor ->", length(nr), "modules,", sum(nr), "rows\n")

qsave(meta, abs(rel("store_meta.qs")))
cat("WROTE", abs(rel("store_meta.qs")), "\n")

## ---- validation: fst reads must equal the source ---------------------------
stopifnot(isTRUE(all.equal(
  as.matrix(read_fst(abs(meta$gp$SEAAD2024$means_path))[, meta$gp$SEAAD2024$celltypes]),
  unname(read_mat(file.path(APP, "www/gene_projection/SEAAD2024/PFC/means/genomewide_means.csv"))),
  check.attributes = FALSE)))
cat("VALIDATION OK — fst store matches source. Done.\n")

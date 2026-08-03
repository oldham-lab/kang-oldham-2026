# Truncated version of project_rand_and_calculate_euclidean:
# Stops after computing projdistpvalindiv and saves it as a .csv file.
# randinds is passed in memory rather than saved/loaded between the two steps.

project_rand_and_calculate_pval <- function(module_output_dir,
                                            filter_under = 3,
                                            do_log = T,
                                            bulk_genes,
                                            save_dir1,
                                            sn_objs,
                                            proj_all,
                                            rand_n = 10000,
                                            seed = 26,
                                            out_prefix = "",
                                            save_randinds = FALSE) {

  logstr <- ifelse(do_log, "_log", "")

  set.seed(seed)

  require(data.table)

  # ── project_random_modules_dt ───────────────────────────────────────────────

  cat("Loading modules and data...\n")
  modules <- fread(file = file.path(module_output_dir, "kme_tables", "topmodposbc_table.csv"))
  modules_list <- split(modules[[2]], modules[[3]])  # Faster than tapply

  modulelengths <- lengths(modules_list)
  these_mods <- which(modulelengths > filter_under)
  modstrimlength <- modulelengths[these_mods]
  sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
  these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]
  
  sn_genes <- rownames(sn_objs$mean[[1]][[1]])

  cat("Generating random modules...\n")
  randmods <- lapply(unique(modstrimlength), function(x) {
    replicate(rand_n, sample(bulk_genes, x), simplify = FALSE)
  })

  cat("Projecting random modules...\n")

  # Ultra-fast projection using matrix operations
  proj_rand_fast <- function(modules, all_m, sn_genes) {
    n_modules <- length(modules)
    n_celltypes <- ncol(all_m)

    # Convert to matrix for speed
    all_m_mat <- as.matrix(all_m)

    # Pre-allocate
    result <- matrix(0, nrow = n_modules, ncol = n_celltypes,
                     dimnames = list(NULL, colnames(all_m)))

    # Vectorized computation
    for (k in seq_along(modules)) {
      idx <- which(sn_genes %in% modules[[k]])
      if (length(idx) > 2) {
        result[k, ] <- colMeans(all_m_mat[idx, , drop = FALSE])
      }
    }

    return(as.data.frame(result))
  }

  randinds <- list()

  for (ct in seq_along(sn_objs$mean)) {
    randinds2 <- list()

    for (split in seq_along(sn_objs$mean[[ct]])) {
      cat(sprintf("  Processing ct=%d, split=%d\n", ct, split))

      all_m <- sn_objs$mean[[ct]][[split]]

      randinds2[[split]] <- future_lapply(
        randmods,
        proj_rand_fast,
        all_m = all_m,
        sn_genes = sn_genes,
        future.globals = structure(TRUE, add = c("all_m", "sn_genes"))
      )
    }

    names(randinds2) <- names(sn_objs$mean[[ct]])
    randinds[[ct]] <- randinds2
  }

  names(randinds) <- names(sn_objs$mean)

  if (save_randinds) {
    save_dir <- file.path(save_dir1, "rand_mod_dists")
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    qsave(randinds, file = file.path(save_dir, paste0(out_prefix, "randinds_native", logstr, ".qs")))
  }

  cat("Done projecting random modules!\n")

  # ── calculate_euclidean_distances ───────────────────────────────────────────

  # randinds passed directly from above (no qread needed)
  randinds <- unlist(randinds, recursive = F)

  # Calculate euclidean distances for random modules (for individual celltypes)
  cat("Calculating differences between random modules for individual CCs...\n")
  diffabs <- lapply(1:(length(randinds) / 2), \(j) {
    mapply(\(x, y) abs(x - y), randinds[[(j * 2 - 1)]], randinds[[j * 2]], SIMPLIFY = F)
  })

  # Calculate euclidean distances for real modules (for individual celltypes)
  cat("Calculating real euclidean distances for indiv cell classes...\n")
  proj_dist_indiv <- lapply(1:(length(proj_all) / 2), \(j) {
    abs(proj_all[[j * 2 - 1]] - proj_all[[j * 2]])
  })

  # Calculate distance p-values for real modules (for individual celltypes)
  cat("Calculating p-values for indiv cell classes...\n")
  projdistpvalindiv <- lapply(seq_along(proj_dist_indiv), \(x) {
    outlist <- list()
    for (b in 1:nrow(proj_dist_indiv[[x]])) {
      ind <- which(unique(modstrimlength) == modstrimlength[b])
      subc <- c()
      for (c in seq_along(proj_dist_indiv[[x]])) {
        subc[c] <- 1 - sum(proj_dist_indiv[[x]][b, c] > diffabs[[x]][[ind]][, c]) / length(diffabs[[x]][[ind]][, c])
      }
      outlist[[b]] <- subc
    }
    outlist <- do.call(rbind, outlist) %>% as.data.frame
    colnames(outlist) <- colnames(proj_dist_indiv[[x]])
    return(outlist)
  })

  # Save projdistpvalindiv dataframe as .csv
  if (!dir.exists(save_dir1))
    dir.create(save_dir1, recursive = T)
  out_df <- projdistpvalindiv[[1]]
  rownames(out_df) <- 1:nrow(out_df)
  write.csv(out_df, file = file.path(save_dir1, paste0(out_prefix, "projdistpvalindiv.csv")))
}

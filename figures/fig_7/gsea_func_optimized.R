# GSEA related functions (optimized)

library(GSEABase)
library(limma)
library(qvalue)
library(ellipse)
library(WGCNA)
library(flashClust)
library(data.table)
library(parallel)

# ── helpers ──────────────────────────────────────────────────────────────────

.clean_genes <- function(genes) {
  genes <- genes[!is.na(genes) & genes != ""]
  unique(toupper(gsub(" ", "", genes)))
}

# phyper-based one-sided Fisher's exact test (equivalent, much faster)
# all_n : total number of background genes (scalar, passed via closure)
.make_fisher <- function(all_n) {
  function(set, mod) {
    mod_n  <- length(mod)
    shared <- length(intersect(mod, set))
    phyper(shared - 1L,
           length(set),
           all_n - length(set),
           mod_n,
           lower.tail = FALSE)
  }
}

# ── GSHG ─────────────────────────────────────────────────────────────────────

GSHG <- function(allModules,
                 allgenes,
                 set_list,   # a list of genesets with names
                 file_desc,
                 save_name,
                 n_cores = detectCores() - 1L) {

  mySetNames <- names(set_list)

  # Pre-process gene sets once
  mySets <- lapply(set_list, function(x) toupper(as.character(unlist(x))))

  # Pre-process background genes once
  allgenes <- .clean_genes(allgenes)
  all_n    <- length(allgenes)

  # Filter sets to those overlapping any module gene, then intersect with background
  AllModGenes <- unique(toupper(unlist(allModules)))
  overlap     <- sapply(mySets, function(s) length(intersect(s, AllModGenes)))
  keep        <- overlap > 0
  mySets      <- lapply(mySets[keep], function(s) intersect(s, allgenes))
  mySetNames  <- mySetNames[keep]
  file_desc   <- file_desc[keep]

  # Drop sets that are empty after background intersection
  non_empty  <- sapply(mySets, length) > 0
  mySets     <- mySets[non_empty]
  mySetNames <- mySetNames[non_empty]
  file_desc  <- file_desc[non_empty]

  fisherTest <- .make_fisher(all_n)

  # Parallelize over modules (outer loop); sets are fast with phyper
  results_list <- mclapply(allModules, function(mod) {
    mod <- unique(toupper(unlist(mod)))
    sapply(mySets, fisherTest, mod = mod)
  }, mc.cores = n_cores)

  GSHGresults <- do.call(cbind, results_list)
  GSHGresults[GSHGresults == 0] <- 1e-300
  colnames(GSHGresults) <- names(allModules)

  datout <- data.frame(
    "SetID"   = mySetNames,
    "SetName" = file_desc,
    "SetSize" = sapply(mySets, length),
    as.data.frame(GSHGresults),
    stringsAsFactors = FALSE
  )

  # fwrite(datout, file = paste0(save_dir, "/", save_name, ".csv"))
  return(datout)
}

# ── BroadGSHG ────────────────────────────────────────────────────────────────

BroadGSHG <- function(allModules,
                      allgenes,
                      allSets,
                      save_name,
                      n_cores = detectCores() - 1L) {

  # Pre-process background genes once
  allgenes <- .clean_genes(allgenes)
  all_n    <- length(allgenes)

  # Filter sets to those overlapping any module gene
  AllModGenes <- unique(toupper(unlist(allModules)))
  overlap     <- sapply(geneIds(allSets), function(s) length(intersect(s, AllModGenes)))
  allSets     <- allSets[overlap > 0]

  # Extract metadata
  setIDs       <- sapply(allSets, setIdentifier)
  collType     <- lapply(allSets, collectionType)
  catType      <- sapply(collType, bcCategory)
  descriptions <- sapply(allSets, description)
  pubmed       <- as.character(sapply(allSets, pubMedIds))
  pubmed[pubmed == ""] <- NA

  # Pre-intersect set genes with background once
  setGenes <- lapply(allSets, function(s) intersect(toupper(geneIds(s)), allgenes))
  setSize  <- sapply(setGenes, length)

  # Drop empty sets after background intersection
  non_empty <- setSize > 0
  allSets   <- allSets[non_empty]
  setGenes  <- setGenes[non_empty]
  setSize   <- setSize[non_empty]
  setIDs    <- setIDs[non_empty]
  descriptions <- descriptions[non_empty]
  pubmed    <- pubmed[non_empty]

  fisherTest <- .make_fisher(all_n)

  # Parallelize over modules
  results_list <- mclapply(allModules, function(mod) {
    mod <- unique(toupper(unlist(mod)))
    sapply(setGenes, fisherTest, mod = mod)
  }, mc.cores = n_cores)

  GSHGresults <- do.call(cbind, results_list)
  GSHGresults[GSHGresults == 0] <- 1e-300
  colnames(GSHGresults) <- names(allModules)

  datout <- data.frame(
    "SetID"   = setIDs,
    "SetName" = names(allSets),
    "SetSize" = setSize,
    as.data.frame(GSHGresults),
    stringsAsFactors = FALSE
  )

  # fwrite(datout, file = paste0(save_dir, "/", save_name, ".csv"))
  return(datout)
}

# ── run_gsea_for_proj ─────────────────────────────────────────────────────────

run_gsea_for_proj_optimized <- function(allModules,
                              set_list  = NULL,
                              file_desc = NULL,
                              broad     = TRUE,
                              n_cores   = detectCores() - 1L) {

  # ── result cache ────────────────────────────────────────────────────────────
  # GSEA here is a pure function of the module set and the geneset sources, but
  # the Broad pass alone parses a 201 MB MSigDB XML and runs across ~1000 modules.
  # fig_6/v6 and fig_8/v4 build IDENTICAL module sets (same kME table, same
  # filter_under = 3), so uncached, each pays that cost in full -- even when the
  # only thing changing is a panel title. Key on the module contents plus the
  # size+mtime of every geneset source that feeds the chosen branch, so edits to
  # any input miss the cache rather than silently returning stale results.
  # Escape hatch: GSEA_CACHE=0 in the environment bypasses read and write.
  .cache_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_cache")
  .p <- list(
    combat  = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"),
    msigdb  = Sys.getenv("MSIGDB_XML", "/home/gugene/code/git/CoPA/data-raw/msigdb_v7.4.xml"),
    projset = system.file("extdata", "proj_sets.qs", package = "CoPA"),
    setdir  = system.file("extdata", "GeneSets_consensus", package = "CoPA")
  )
  .use_cache  <- !identical(Sys.getenv("GSEA_CACHE"), "0")
  .cache_file <- NULL
  if (.use_cache) {
    .stamp <- function(f) {
      i <- file.info(f)
      paste(f, i$size, format(i$mtime, "%Y%m%d%H%M%S"), sep = "|")
    }
    # Only the sources the chosen branch actually reads: set_list/proj_sets are
    # unused when broad = TRUE, so including them would fragment the cache.
    .src <- .p$combat
    if (broad) {
      .src <- c(.src, .p$msigdb)
    } else {
      .src <- c(.src, .p$projset,
                if (is.null(set_list)) sort(list.files(.p$setdir, full.names = TRUE)) else character(0))
    }
    .key <- digest::digest(list(
      modules  = lapply(allModules[order(names(allModules))], function(g) sort(as.character(g))),
      modnames = sort(names(allModules)),
      broad    = broad,
      setlist  = if (broad || is.null(set_list)) NULL
                 else lapply(set_list, function(x) sort(as.character(unlist(x)))),
      filedesc = if (broad) NULL else file_desc,
      sources  = vapply(.src, .stamp, character(1), USE.NAMES = FALSE)
    ), algo = "xxhash64")
    .cache_file <- file.path(.cache_dir, paste0("gsea_",
                             if (broad) "broad" else "userInput", "_", .key, ".qs"))
    if (file.exists(.cache_file)) {
      cat("GSEA cache hit:", basename(.cache_file), "\n")
      return(qs::qread(.cache_file))
    }
  }

  allgenes <- fread(
    file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"),
    data.table = FALSE
  )[, 2]

  if (is.null(set_list)) {
    set_dir    <- system.file("extdata", "GeneSets_consensus", package = "CoPA")
    gsea_legend <- read.csv(paste0(set_dir, "/ConsensusGeneSets.csv"))
    files      <- list.files(set_dir, full.names = TRUE)[grep("MOSET", list.files(set_dir))]
    file_names <- list.files(set_dir)[grep("MOSET", list.files(set_dir))]

    set_list  <- vector("list", length(files))
    file_desc <- character(length(files))
    for (l in seq_along(files)) {
      set_list[[l]]      <- read.csv(files[l], header = FALSE)
      names(set_list)[l] <- gsub(".csv", "", file_names[l])
      file_desc[l]       <- as.character(
        gsea_legend[gsea_legend[, 1] == names(set_list)[[l]], 2]
      )
    }
  }

  set_list_sn <- qread(system.file("extdata", "proj_sets.qs", package = "CoPA"))
  file_desc_sn <- names(set_list_sn)

  .res <- if (broad) {
    broadSets <- getBroadSets(.p$msigdb)
    cat("Running GSEA using Broad genesets...\n")
    BroadGSHG(
      allModules = allModules,
      allgenes   = allgenes,
      allSets    = broadSets,
      save_name  = "gsea_Broad",
      n_cores    = n_cores
    )
  } else {
    cat("Running GSEA using user input genesets...\n")
    GSHG(
      allModules = allModules,
      allgenes   = allgenes,
      set_list   = c(set_list, set_list_sn),
      file_desc  = c(file_desc, file_desc_sn),
      save_name  = "gsea_userInput",
      n_cores    = n_cores
    )
  }

  if (.use_cache) {
    dir.create(.cache_dir, showWarnings = FALSE, recursive = TRUE)
    qs::qsave(.res, .cache_file)
    cat("GSEA cached:", basename(.cache_file), "\n")
  }
  .res
}

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

  if (broad) {
    broadSets <- getBroadSets(Sys.getenv("MSIGDB_XML", "/home/gugene/code/git/COPA/data/msigdb_v7.4.xml"))
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
}

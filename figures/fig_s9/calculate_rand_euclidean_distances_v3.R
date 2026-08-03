# Modification of COPA/calculate_euclidean_distances.R
# create a "network" of random modules, with mod lengths matching real modules,
# then calculate p-values (after taking rand mods out of rand distribution)
# v3: proj_all and randinds2 passed directly (loading moved to calling script)

calculate_euclidean_distances <- function(real, # T or F
                                          proj_all = NULL,
                                          randinds2,
                                          module_output_dir,
                                          filter_under = 3,
                                          seed = 23,
                                          save_dir1,
                                          file_index = "",
                                          save_direction = F,
                                          save_pval = F){

  if(!dir.exists(save_dir1))
    dir.create(save_dir1, recursive = T)

  # Load log native indices and other necessary objects
  modulelengths <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv")) |>
      (\(.) tapply(.[,2], .[,3], list))() |>
      lapply(length) |> unlist()
  these_mods <- which(modulelengths>filter_under)
  #sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
  #these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)] # need to edit upstream scripts to account for this
  modstrimlength <- modulelengths[these_mods]

  if(!real){
    # Create a "network" of random modules
    set.seed(seed)
    randindlist <- lapply(seq_along(modstrimlength), \(x) sample(1:10000, 1))
    proj_all <- lapply(1:2, \(n){
      mapply(\(x, y){
        which.rand <- which(unique(modstrimlength) == modstrimlength[x])
        return(randinds2[[n]][[which.rand]][y, ])
      }, seq_along(modstrimlength), randindlist, SIMPLIFY = F) |>
        do.call(what = "rbind")
    })
  }

  # Calculate euclidean distances for random modules (for individual celltypes)
  # (subclass only)
  cat("Calculating differences between random modules for individual CCs...\n")
  diffabs <- lapply(1:(length(randinds2)/2), \(j) {
    bla <- mapply(\(x,y){
      abs(x-y)
    },randinds2[[(j*2-1)]], randinds2[[j*2]], SIMPLIFY=F)
  })

  if(!real){
    # Remove modules used to create random "network" from list
    # (subclass only)
    for(i in 1:length(modstrimlength)){
      which.rand <- which(unique(modstrimlength) == modstrimlength[i])
      diffabs[[1]][[which.rand]] <- diffabs[[1]][[which.rand]][-randindlist[[i]], ]
    }
  }

  # Calculate euclidean distances for individual celltypes
  cat("Calculating real euclidean distances for indiv cell classes...\n")
  proj_dist_indiv <- lapply(1:(length(proj_all)/2), \(j){
    abs(proj_all[[j*2-1]] - proj_all[[j*2]])
  })

  # Save direction matrix
  if(save_direction){
    proj_dist_dir <- lapply(1:(length(proj_all)/2), \(j){
      proj_all[[j*2-1]] - proj_all[[j*2]]
    })
    fwrite(proj_dist_dir[[1]], file = paste0(save_dir1, "/mod_direction", file_index, ".csv"))
  }


  # Calculate distance p-values for individual celltypes
  cat("Calculating p-values for indiv cell classes...\n")
  projdistpvalindiv <- lapply(seq_along(proj_dist_indiv), \(x){
    outlist <- list()
    for(b in 1:nrow(proj_dist_indiv[[x]])){
      ind <- which(unique(modstrimlength) == modstrimlength[b])
      subc <- c()
      for(c in seq_along(proj_dist_indiv[[x]])){
        subc[c] <- 1 - sum(proj_dist_indiv[[x]][b,c] > diffabs[[x]][[ind]][,c]) / length(diffabs[[x]][[ind]][,c])
      }
      outlist[[b]] <- subc
    }
    outlist <- do.call(rbind, outlist) %>% as.data.frame
    colnames(outlist) <- colnames(proj_dist_indiv[[x]])
    return(outlist)
  })

  # Save (subclass only)
  if(save_pval){
    if(real){
      fwrite(projdistpvalindiv[[1]], file = paste0(save_dir1, "/real_pvals.csv"))
    } else {
      fwrite(projdistpvalindiv[[1]], file = paste0(save_dir1, "/rand_pvals", file_index, ".csv"))
    }
  } else {
    return(projdistpvalindiv[[1]])
  }
}

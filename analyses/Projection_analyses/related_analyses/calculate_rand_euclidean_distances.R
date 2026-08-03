# Modification of COPA/calculate_euclidean_distances.R
# create a "network" of random modules, with mod lengths matching real modules, 
# then calculate p-values (after taking rand mods out of rand distribution)

calculate_euclidean_distances <- function(real, # T or F
                                          sn_summary_object_path,
                                          module_output_dir,
                                          filter_under = 3,
                                          logstr = "_log",
                                          seed = 23,
                                          save_dir1,
                                          file_index = "",
                                          randinds2 = NULL,
                                          save_direction = F){

  if(!dir.exists(save_dir1))
    dir.create(save_dir1, recursive = T)
  
  # Load log native indices and other necessary objects
  modulelengths <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv")) |>
      (\(.) tapply(.[,2], .[,3], list))() |>
      lapply(length) |> unlist()
  these_mods <- which(modulelengths>filter_under)
  modstrimlength <- modulelengths[these_mods]

  if(is.null(randinds2)){
    randinds2 <- qread(file.path(sn_summary_object_path,"rand_mod_dists","randinds_native_log.qs")) |>
        unlist(recursive=F)
    randinds2 <- randinds2[1:2] # subclass only
  }

  if(real){
    sn_objs <- qread(file = file.path(sn_summary_object_path, "sn_summary_tables",paste0("sn_summary_objects",logstr,".qs")))
    snnames <- mapply(\(ct, ct_name) paste0(ct_name, "_", names(ct)), sn_objs$mean, names(sn_objs$mean), SIMPLIFY=F) 
    proj_all <- lapply(unlist(snnames)[1:2], \(x){ # subclass only
      index_dir <- list.files(file.path(sn_summary_object_path, "sn_proj_indices","log_native"),full.names=T)
      files1 <- index_dir[grep("indices_over",index_dir)]
      fread(data.table=F,file=files1[grep(x,files1)])[these_mods,]
    })
  } else {
    # Create a "network" of random modules
    set.seed(seed)
    randindlist <- lapply(seq_along(modstrimlength), \(x) sample(1:10000, 1))
    proj_all <- lapply(1:2, \(i){
      mapply(\(x, y){
        which.rand <- which(unique(modstrimlength) == modstrimlength[x])
        return(randinds2[[i]][[which.rand]][y, ])
      }, seq_along(modstrimlength), randindlist, SIMPLIFY = F) |> 
        do.call(what = "rbind") %>%
        mutate(module = 1:nrow(.))
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
    rem <- which(colnames(proj_all[[j*2-1]])=="module")
    return(abs(proj_all[[j*2-1]][,-rem] - proj_all[[j*2]][,-rem]))
  })

  # Save direction matrix
  if(save_direction){
    proj_dist_dir <- lapply(1:(length(proj_all)/2), \(j){
      rem <- which(colnames(proj_all[[j*2-1]])=="module")
      proj_all[[j*2-1]][,-rem] - proj_all[[j*2]][,-rem]
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
  if(real){
    fwrite(projdistpvalindiv[[1]], file = paste0(save_dir1, "/real_pvals.csv"))
  } else {
    fwrite(projdistpvalindiv[[1]], file = paste0(save_dir1, "/rand_pvals", file_index, ".csv"))
  }
}




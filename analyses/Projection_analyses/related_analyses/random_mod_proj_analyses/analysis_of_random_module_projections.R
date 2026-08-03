library(data.table)
library(qs)
library(tidyverse)
# 1. create some random mods
# 2. project on sn datasets
# 3. evaluate random 

######################
# Calculate random modules (run once)
######################

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)
mod_seed_lengths <- lapply(mod_seed,length) |> unlist()

# Make 100000 rand modules
# - sample from mod seed length dist
set.seed(23)
randmods <- lapply(1:100000, \(x){
  sample(bulk_expr$Gene, size = sample(mod_seed_lengths, 1))
})
qsave(randmods, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/randmods.qs"))

# Calculate inds for each of 4 datasets
# Simplified index function
calculate_proj_indices2 <- function(norm_by_mean=F,
                                   norm_by_max=F,
                                   save_dir,
                                   modules,
                                   sn_summary_object_path,
                                   logstr = "_log"
){

  if(!dir.exists(save_dir))
    dir.create(save_dir, recursive = T)

  sn_objs <- qread(file = file.path(sn_summary_object_path, "sn_summary_tables",paste0("sn_summary_objects",logstr,".qs")))

  for(ct in seq_along(sn_objs$mean)){
    for(split in seq_along(sn_objs$mean[[ct]])){
      all_m <- sn_objs$mean[[ct]][[split]]
      all_var <- sn_objs$var[[ct]][[split]]
      var_na <- apply(all_var,2,function(x) sum(is.na(x)))
      all_n <- sn_objs$n[[ct]][[split]]
      sn_name <- paste0(names(sn_objs$mean)[ct], "_", names(sn_objs$mean[[ct]])[split])
      all_means <- sn_objs$all_means[[ct]]
      
      cat(paste0("Calculating projection indices for ", sn_name,"...\n"))
      
      # for each module, find mean across all mod genes
      proj_ind <- list()
     # proj_se <- list()
      for(k in 1:length(modules)){ # for each module
        # subset to module genes
        boolvec <- rownames(all_m) %in% modules[[k]]
        
        if(sum(boolvec) > 1){ # if there is more than one gene in module
          m_trim <- all_m[rownames(all_m) %in% modules[[k]], ]
          var_trim <- all_var[rownames(all_var) %in% modules[[k]], ]
          cc_mean_out <- list()
          #cc_se_out <- list()
          for(l in 1:ncol(m_trim)){ #for all cell classes
            if(var_na[l]==0){ # if no NAs
              # average means/var across mod genes
              decomp_df <- utilities::sample.decomp(n = rep(all_n[l], nrow(m_trim)),
                                                    sample.mean = m_trim[,l],
                                                    sample.var = var_trim[,l])
              pool <- nrow(decomp_df)
              cc_mean_out[[l]] <- decomp_df[pool,2]
              #cc_se_out[[l]] <- sqrt(decomp_df[pool,3]/decomp_df[pool,1])
              if(norm_by_mean){
                cc_mean_out[[l]] <- cc_mean_out[[l]]/all_means[l]
               # cc_se_out[[l]] <- cc_se_out[[l]]/all_means[l]
              }
            } else {
              cc_mean_out[[l]] <- mean(m_trim[,l])
            # cc_se_out[[l]] <- 0
           }
          } # end of for l
          
          if(norm_by_max){
            cc_mean_out <- lapply(cc_mean_out, function(x) x/max(unlist(cc_mean_out), na.rm=T))
          #  cc_se_out <- lapply(cc_se_out, function(x) x/max(unlist(cc_mean_out), na.rm=T))
          }
          
          names(cc_mean_out) <- colnames(all_m)
        #  names(cc_se_out) <- colnames(all_m)

          proj_ind[[k]] <- unlist(cc_mean_out)
         # proj_se[[k]] <- unlist(cc_se_out)
        } else if(sum(boolvec)==1){
          m_trim <- all_m[rownames(all_m) %in% modules[[k]], ]
          proj_ind[[k]] <- m_trim
        #  proj_se[[k]] <- rep(0, ncol(all_m))
         # names(proj_se[[k]]) <- names(m_trim)
        } else {
          proj_ind[[k]] <- rep(0, ncol(all_m))
          #proj_se[[k]] <- rep(0, ncol(all_m))
          #names(proj_se[[k]]) <- names(m_trim)
        }
      } # for k
      all_df <- as.data.frame(do.call(rbind, proj_ind))
      colnames(all_df)[1:ncol(all_m)] <- colnames(all_m) 
      all_df$module <- names(modules)
      #all_df_se <- as.data.frame(do.call(rbind, proj_se))
      #colnames(all_df_se)[1:ncol(all_m)] <- colnames(all_m)
      #all_df_se$module <- names(modules)
      
      write.csv(all_df, file = file.path(save_dir, paste0("indices_over_all_datasets_",sn_name,".csv")), row.names=F)
      #write.csv(all_df_se, file = file.path(save_dir, paste0("indices_se_over_all_datasets_",sn_name,".csv")), row.names=F)
    }
  }
}

## Calculate indices for random modules:
# Output dir
save_dir1 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/")

# SEA
calculate_proj_indices2(save_dir = file.path(save_dir1, "sea"),
                        modules = randmods,
                        sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
                        )

# Lein
calculate_proj_indices2(save_dir = file.path(save_dir1, "jorstad"),
                        modules = randmods,
                        sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinMTG/")
                        )

# MIT
calculate_proj_indices2(save_dir = file.path(save_dir1, "mit"),
                        modules = randmods,
                        sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/")
                        )

# Morabito
calculate_proj_indices2(save_dir = file.path(save_dir1, "morabito"),
                        modules = randmods,
                        sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/Morabito_ABIanno/")
                        )

## Load indices, align celltypes, calculate upper tri of index pairwise cors
base_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/")
index_dirs <- c("log_native", "log_REI")
index_list <- list()
index_save_names <- c("native_log", "REI")
index_xaxis_names <- c("Mean expression (log UMI counts)", "Relative expression index")
#for(d in 1:2){
d=1
  lein <- fread(data.table = F, file = paste0(base_dir, "/jorstad/indices_over_all_datasets_Subclass_con.csv")) 

  mit <- fread(data.table = F, file = paste0(base_dir, "/mit/indices_over_all_datasets_Subclass_Control.csv")) 

  SEAcon <- fread(data.table = F, file = paste0(base_dir, "/sea/indices_over_all_datasets_Subclass_Control.csv")) 

  mora <- fread(data.table = F, file = paste0(base_dir, "/morabito/indices_over_all_datasets_Subclass_Control.csv")) 
  # Change celltypes to match

  colnames(lein)[15] <- "Micro"
  colnames(mit) <- gsub("Exc ", "", colnames(mit))
  colnames(mit) <- gsub("Inh ", "", colnames(mit))
  colnames(mit) <- gsub("-", "/", colnames(mit))
  colnames(mit)[c(2, 3, 4, 6, 11, 12, 16, 21, 25, 27)] <- c("Micro", "Vip", "Astro", "Oligo", "Pvalb", "Endo", "Sst", "Lamp5", "Pax6", "L6 IT Car3")
  colnames(mora)[c(1, 3, 14, 15)] <- c("Astro", "Endo", "Micro", "Oligo")
  colnames(SEAcon)[c(1, 3, 15, 16)] <- c("Astro", "Endo", "Micro", "Oligo")

  # Find intersection of celltypes and create euler diagram
  commonct <- Reduce(intersect, list(colnames(lein), colnames(mit), 
                   # colnames(scope), 
                    colnames(SEAcon), colnames(mora)))
  allcts <- unique(c(colnames(lein), colnames(mit), 
                    #colnames(scope), 
                    colnames(SEAcon), colnames(mora)))

  # Hardcode order of cts
  #  "excitatory neurons, inhibitory neurons, astrocytes, oligodendrocytes, OPCs, 
  #  microglia, vascular cells (including fibroblasts), T cells."
  # (remove EC - entorhinal cortex-specific celltype from Liu 2025)
  allcts <- c("L2/3 IT", "L3/4 IT", "L3/5 IT", "L4 IT", "L4/5 IT/1", "L4/5 IT/2", "L5 ET", "L5 IT",  "L5/6 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro", "Endo", "Fib", "Per", "SMC", "VLMC", "T")

  # Fill out celltypes
  lein[setdiff(allcts, colnames(lein))] <- NA
  lein <- lein[, match(allcts, colnames(lein))]

  mit[setdiff(allcts, colnames(mit))] <- NA
  mit <- mit[, match(allcts, colnames(mit))]

  SEAcon[setdiff(allcts, colnames(SEAcon))] <- NA
  SEAcon <- SEAcon[, match(allcts, colnames(SEAcon))]

  mora[setdiff(allcts, colnames(mora))] <- NA
  mora <- mora[, match(allcts, colnames(mora))]

  # Plot individual cor heatmaps for each module
  corlistind <- list()
  for(i in 1:nrow(lein)){
    tempdf <- data.frame("Jorstad 2023" = unlist(lein[i, ]), "Liu 2025" = unlist(mit[i, ]), 
                         #"brainSCOPE_con" = unlist(scope[i, ]), 
                         "Gabitto 2024" = unlist(SEAcon[i, ]), "Morabito 2021" = unlist(mora[i, ]))
    tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
    #colnames(tempcor) <- c("Jorstad 2023",  "Liu 2025", "Gabitto 2024", "Morabito 2021")
    #rownames(tempcor) <- c("Jorstad 2023",  "Liu 2025", "Gabitto 2024", "Morabito 2021")
    # testvec <- apply(tempcor, 2, \(x) sum(is.na(x)))
    # if(any(testvec > 0)){
    #   rem <- which(testvec > 0)
    #   tempcor <- tempcor[-rem, -rem]
    # } 
    corlistind[[i]] <- tempcor[upper.tri(tempcor)]
  }
#}

qsave(corlistind, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/randind_nativelog.qs"))

## Load real indices, align, calculate upper tri of index pairwise cors
base_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")
index_dirs <- c("log_native", "log_REI")
index_list <- list()
index_save_names <- c("native_log", "REI")
index_xaxis_names <- c("Mean expression (log UMI counts)", "Relative expression index")
#for(d in 1:2){
d=1
  lein <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Cell_Type_1.csv")) |>
    dplyr::select(!module)
  lein_se <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_se_Cell_Type.csv")) 
  lein_se <- lein_se[, match(colnames(lein), colnames(lein_se))]

  mit <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
    dplyr::select(!module)
  mit_se <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_se_RNA.SubclassCon.csv")) 
  mit_se <- mit_se[, match(colnames(mit), colnames(mit_se))]

  SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
    dplyr::select(!module)
  SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassCon.csv")) 
  SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]

  mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
    dplyr::select(!module)
  mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassCon.csv")) 
  mora_se <- mora_se[, match(colnames(mora), colnames(mora_se))]

  # Change celltypes to match

  colnames(lein)[15] <- "Micro"
  colnames(mit) <- gsub("Exc ", "", colnames(mit))
  colnames(mit) <- gsub("Inh ", "", colnames(mit))
  colnames(mit) <- gsub("-", "/", colnames(mit))
  colnames(mit)[c(2, 3, 4, 6, 11, 12, 16, 21, 25, 27)] <- c("Micro", "Vip", "Astro", "Oligo", "Pvalb", "Endo", "Sst", "Lamp5", "Pax6", "L6 IT Car3")
  colnames(mora)[c(1, 3, 14, 15)] <- c("Astro", "Endo", "Micro", "Oligo")
  colnames(SEAcon)[c(1, 3, 15, 16)] <- c("Astro", "Endo", "Micro", "Oligo")

  colnames(lein_se) <- colnames(lein)
  colnames(mit_se) <- colnames(mit)
  colnames(SEAcon_se) <- colnames(SEAcon)
  colnames(mora_se) <- colnames(mora)

  # Find intersection of celltypes and create euler diagram
  commonct <- Reduce(intersect, list(colnames(lein), colnames(mit), 
                   # colnames(scope), 
                    colnames(SEAcon), colnames(mora)))
  allcts <- unique(c(colnames(lein), colnames(mit), 
                    #colnames(scope), 
                    colnames(SEAcon), colnames(mora)))

  # Hardcode order of cts
  #  "excitatory neurons, inhibitory neurons, astrocytes, oligodendrocytes, OPCs, 
  #  microglia, vascular cells (including fibroblasts), T cells."
  # (remove EC - entorhinal cortex-specific celltype from Liu 2025)
  allcts <- c("L2/3 IT", "L3/4 IT", "L3/5 IT", "L4 IT", "L4/5 IT/1", "L4/5 IT/2", "L5 ET", "L5 IT",  "L5/6 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro", "Endo", "Fib", "Per", "SMC", "VLMC", "T")

  # Fill out celltypes
  lein[setdiff(allcts, colnames(lein))] <- NA
  lein <- lein[, match(allcts, colnames(lein))]
  lein_se[setdiff(allcts, colnames(lein_se))] <- NA
  lein_se <- lein_se[, match(allcts, colnames(lein_se))]

  mit[setdiff(allcts, colnames(mit))] <- NA
  mit <- mit[, match(allcts, colnames(mit))]
  mit_se[setdiff(allcts, colnames(mit_se))] <- NA
  mit_se <- mit_se[, match(allcts, colnames(mit_se))]

  SEAcon[setdiff(allcts, colnames(SEAcon))] <- NA
  SEAcon <- SEAcon[, match(allcts, colnames(SEAcon))]
  SEAcon_se[setdiff(allcts, colnames(SEAcon_se))] <- NA
  SEAcon_se <- SEAcon_se[, match(allcts, colnames(SEAcon_se))]

  mora[setdiff(allcts, colnames(mora))] <- NA
  mora <- mora[, match(allcts, colnames(mora))]
  mora_se[setdiff(allcts, colnames(mora_se))] <- NA
  mora_se <- mora_se[, match(allcts, colnames(mora_se))]

  # Plot individual cor heatmaps for each module
  corlistind <- list()
  for(i in 1:nrow(lein)){
    tempdf <- data.frame("Jorstad 2023" = unlist(lein[i, ]), "Liu 2025" = unlist(mit[i, ]), 
                         #"brainSCOPE_con" = unlist(scope[i, ]), 
                         "Gabitto 2024" = unlist(SEAcon[i, ]), "Morabito 2021" = unlist(mora[i, ]))
    tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
    # colnames(tempcor) <- c("Jorstad 2023",  "Liu 2025", "Gabitto 2024", "Morabito 2021")
    # rownames(tempcor) <- c("Jorstad 2023",  "Liu 2025", "Gabitto 2024", "Morabito 2021")
    # testvec <- apply(tempcor, 2, \(x) sum(is.na(x)))
    # if(any(testvec > 0)){
    #   rem <- which(testvec > 0)
    #   tempcor <- tempcor[-rem, -rem]
    # } 

    corlistind[[i]] <- tempcor[upper.tri(tempcor)]
  }
#}
qsave(corlistind, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/realind_nativelog.qs"))

###############
#New session
###############

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)
mod_seed_lengths <- lapply(mod_seed,length) |> unlist()

# Load random mods
randmods <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/randmods.qs"))

# Load real and random index pairwise cors (upper.tri)
corlistind <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/realind_nativelog.qs"))
corlistrand <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/randind_nativelog.qs"))


# Plot distribution of real cors vs rand
plotdf <- list(data.frame("cor" = lapply(corlistind, mean)|> unlist(),
                          "type" = "Real modules"),
               data.frame("cor" = lapply(corlistrand, mean) |> unlist(),
                          "type" = "Random modules")) |>
  do.call("what" = "rbind")
plotdf |> group_by(type) |> summarise(meancor = mean(cor, na.rm = T),
                                      medcor = median(cor, na.rm = T))
#   type           meancor
#   <chr>            <dbl>
# 1 Random modules   0.940
# 2 Real modules     0.932

wilcox.test(plotdf$cor[plotdf$type == "Real modules"], plotdf$cor[plotdf$type=="Random modules"], alternative = "greater")
# p-value = 0.001157

p <- ggplot(plotdf, aes(x = cor, color = type)) + 
  theme_classic() +
  geom_density() +
  labs(x = "Mean pairwise cors of projection indices", 
       y = "Density",
       title = "Mean pairwise cors of CoPA projections",
       subtitle = "Projection onto 4 SN-RNAseq datasets\n100,000 random mods vs 1158 real mods\nDistribution of lengths of rand mods matches real\np = 0.0012") +
  theme(legend.title = element_blank(),
        plot.subtitle = element_text(size = 8),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.box.background = element_rect(color = "black", linewidth = 0.5)) +
  geom_vline(xintercept = 0.947, color = "#F8766D", linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0.950, color = "#00BFC4", linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 0.89, y = 20, color = "#F8766D", label = "Median:\n0.947") +
  annotate("text", x = 1, y = 18, color = "#00BFC4", label = "Median:\n0.950") 

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/real_vs_rand_proj_cor.png"), height = 4, width = 6)  
 
# Plot module cors vs projection cors
bulkcor <- cor(t(bulk_expr[,-c(1:2)])) # takes a few minutes
colnames(bulkcor) <- bulk_expr$Gene
rownames(bulkcor) <- bulk_expr$Gene

# Mean pairwise cors of real mods
rbc <- lapply(mod_seed, \(x){
  t <- colnames(bulkcor) %in% x
  temp <- bulkcor[t, t]
  return(temp[upper.tri(temp)])
}) |> 
  lapply(mean) |> 
  unlist()

rabc <- lapply(randmods, \(x){
  t <- colnames(bulkcor) %in% x
  temp <- bulkcor[t, t]
  return(temp[upper.tri(temp)])
}) |>
  lapply(mean) |>
  unlist()

plotdf2 <- plotdf |>
  mutate(bulkcor = c(rbc, rabc))

p2 <- ggplot(plotdf2, aes(y = cor, x = bulkcor, color = type)) + 
  theme_classic() +
  geom_point(alpha = 0.01) +
  geom_smooth() +
  labs(y = "Mean pairwise cors of projection indices", 
       x = "Mean pairwise cors of module seed genes",
       title = "Comparison of mean pairwise cors of seed genes vs projections",
       subtitle = "Projection onto 4 datasets, 10000 random mods vs 1158 real mods\nDistribution of lengths of rand mods matches real") +
  theme(legend.title = element_blank(),
        plot.subtitle = element_text(size = 8),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.box.background = element_rect(color = "black", linewidth = 0.5)) +
  annotate("text", x = 0.2, y = 0.8, label = "r = 0.20", color = "#F8766D") +
  annotate("text", x = 0.7, y = 0.8, label = "r = 0.35", color = "#00BFC4")
ggsave(p2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/projcor_vs_modcor.png"), height = 4, width = 6)  
# Calculate cors and p-vals
plotdf2 |> group_by(type) |>
  summarise(corsum = cor(cor, bulkcor, use='pairwise.complete.obs'),
            corpval = cor.test(cor, bulkcor)$p.value) 
#   type           corsum  corpval
#   <chr>           <dbl>    <dbl>
# 1 Random modules  0.196 0       
# 2 Real modules    0.346 9.64e-34
cor(plotdf2$cor, plotdf2$bulkcor, use = 'pairwise.complete.obs')
# 0.0607

p3 <- cowplot::plot_grid(p, p2, nrow = 1)
ggsave(p3, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/analysis_of_projections_of_randomMods.png"), height = 4, width = 12)

# Calc % VE

expr_t <- t(bulk_expr[,3:ncol(bulk_expr)])
colnames(expr_t) <- bulk_expr[,2]
expr_z <- apply(expr_t, 2, function(x) (x-mean(x))/sd(x))
r2_vec_bc_list <- list()
for(j in seq_along(mod_seed)){
  expr_plot <- as.data.frame(expr_z[ ,colnames(expr_z) %in% mod_seed[[j]]])
  pc1_bc <-  prcomp(expr_plot[apply(expr_plot,2,var)>0], scale=T)$x[,1]
  r2_vec_bc <- apply(expr_plot[apply(expr_plot,2,var)>0],2,function(x) summary(lm(x~pc1_bc))$r.squared)
  r2_vec_bc_list[[j]] <- r2_vec_bc
  cat(j, " ")
}
vevec <- lapply(r2_vec_bc_list, mean) |> unlist()

# calc %VE for rand mods
r2_vec_rand_list <- list()
for(j in seq_along(randmods)){
  expr_plot <- as.data.frame(expr_z[ ,colnames(expr_z) %in% randmods[[j]]])
  pc1_bc <-  prcomp(expr_plot[apply(expr_plot,2,var)>0], scale=T)$x[,1]
  r2_vec_bc <- apply(expr_plot[apply(expr_plot,2,var)>0],2,function(x) summary(lm(x~pc1_bc))$r.squared)
  r2_vec_rand_list[[j]] <- r2_vec_bc
  cat(j, " ")
}
vevec_rand <- lapply(r2_vec_rand_list, mean) |> unlist()

cor(plotdf |> filter(type=="Real modules") |> pull(cor), vevec, use = 'pairwise.complete.obs')
# 0.34

plotdf3 <- plotdf2 |>
  mutate(ve = c(vevec, vevec_rand))

########################################
# What if we produced these analysis after excluding modules that are not significant in at least 2 datasets?
########################################

sig_7 <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/no_of_sig_bulk_datasets_per_mod.csv"))
pcutfdr <- plotdf2[-which(sig_7[,2] < 2), ] 
## Plot after excluding mods not sig in at least 2 datasets (FDR cutoff)
# Calculate cors first:
pcutfdr |> group_by(type) |>
  summarise(corsum = cor(cor, bulkcor, use='pairwise.complete.obs'),
            corpval = cor.test(cor, bulkcor)$p.value) 
# 1 Random modules  0.196 0         
# 2 Real modules    0.191 0.00000267

pcut2 <- ggplot(pcutfdr, aes(y = cor, x = bulkcor, color = type)) + 
  theme_classic() +
  geom_point(alpha = 0.01) +
  geom_smooth() +
  labs(y = "Mean pairwise cors of projection indices", 
       x = "Mean pairwise cors of module seed genes",
       title = "Comparison of mean pairwise cors of seed genes vs projections",
       subtitle = "Projection onto 4 datasets, 10000 random mods vs 807 real mods\nDistribution of lengths of rand mods matches real") +
  theme(legend.title = element_blank(),
        plot.subtitle = element_text(size = 8),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.box.background = element_rect(color = "black", linewidth = 0.5)) +
  annotate("text", x = 0.2, y = 0.8, label = "r = 0.196", color = "#F8766D") +
  annotate("text", x = 0.7, y = 0.8, label = "r = 0.335", color = "#00BFC4")
ggsave(pcut2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/projcor_vs_modcor_filteredforModssignificantIn2OrMoredatasetsFDR.png"), height = 4, width = 6)  


pcutfdr |> group_by(type) |> summarise(meancor = mean(cor, na.rm = T),
                                       medcor = median(cor, na.rm = T))
#   type           meancor
#   <chr>            <dbl>
# 1 Random modules   0.947
# 2 Real modules     0.952

wilcox.test(pcutfdr$cor[pcutfdr$type == "Real modules"], pcutfdr$cor[pcutfdr$type=="Random modules"], alternative = "greater")
# p-value < 2.2e-16

pcut <- ggplot(pcutfdr, aes(x = cor, color = type)) + 
  theme_classic() +
  geom_density() +
  labs(x = "Mean pairwise cors of projection indices", 
       y = "Density",
       title = "Mean pairwise cors of CoPA projections\n(filtered for mods significant in >= 2 datasets)",
       subtitle = "Projection onto 4 SN-RNAseq datasets\n100,000 random mods vs 807 real mods\nDistribution of lengths of rand mods matches real\np < 5.7e-13") +
  theme(legend.title = element_blank(),
        plot.subtitle = element_text(size = 8),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.box.background = element_rect(color = "black", linewidth = 0.5)) +
  geom_vline(xintercept = 0.947, color = "#F8766D", linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0.952, color = "#00BFC4", linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 0.89, y = 20, color = "#F8766D", label = "Median:\n0.947") +
  annotate("text", x = 1, y = 20, color = "#00BFC4", label = "Median:\n0.952") 

ggsave(pcut, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/real_vs_rand_proj_cor_filteredforModssignificantIn2OrMoredatasetsFDR.png"), height = 4, width = 6)  

pall <- cowplot::plot_grid(p, pcut, nrow = 1) |> ggplotify::as.ggplot()
ggsave(pall, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/real_vs_rand_proj_cor_unfilteredVsfiltered.png"), height = 4, width = 12)

library(tidyverse)
library(data.table)
library(qs)
library(ggpubr)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/calculate_rand_euclidean_distances.R"))

# conVAll
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome//"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITMTG/"),
                  caption1 = "Dataset used: Liu et al. 2025 (MIT_Multiome_Multiregion) control vs all AD samples",
                  splits = c("all AD", "con"),
                  calc_rand_mods = T,
                  plot_real_vs_rand = T
                  )
# # conVEarly
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly//"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITMTG_conVEarly/"),
                  caption1 = "Dataset used: Liu et al. 2025 (MIT_Multiome_Multiregion) control vs early AD samples",
                  splits = c("con", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F
                  )

# # earlyVLate
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITMTG_earlyVLate/"),
                  caption1 = "Dataset used: Liu et al. 2025 (MIT_Multiome_Multiregion) early vs late AD samples",
                  splits = c("late AD", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F
                  )

# Function:
rand_mod_analysis <- function(sn_summary_object_path,
                              module_output_dir,
                              save_dir1,
                              caption1,
                              splits,
                              calc_rand_mods,
                              plot_real_vs_rand){

  ###############
  # i) the % of all modules (real or random) that show significant DE in any subclass,  
  ###############                              
  # Calculate dCoPA p-values for real mods

  calculate_euclidean_distances(real = T,
                                sn_summary_object_path = sn_summary_object_path,
                                module_output_dir = module_output_dir,
                                save_dir1 = save_dir1,
                                save_direction = T)
  # qvalue function:
  qval_redund <- \(allPvalues){
    out <- tryCatch({
      return(WGCNA::qvalue(allPvalues))
      }, error = \(cond) {
        message("Adjusted max lambda to max p-value in q-value calculation")
        return(WGCNA::qvalue(allPvalues,lambda = seq(0.05,max(allPvalues), 0.05)))
      },
      warning = \(cond) {
        message(cond)
      }
    )
  }

  # Load real dCoPA p-values for SEAAD2024 (con vs All)
  pvals <- fread(data.table = F, paste0(save_dir1, "real_pvals.csv"))

  if(calc_rand_mods){
    # Calculate dCoPA p-values for random mods
    randinds2 <- qread(file.path(sn_summary_object_path,"rand_mod_dists","randinds_native_log.qs")) |>
        unlist(recursive=F) %>%
        `[`(1:2)

    for(i in 1:1000){ # 1000 iterations takes about an hour
      calculate_euclidean_distances(real = F,
                                    sn_summary_object_path = sn_summary_object_path,
                                    module_output_dir = module_output_dir,
                                    seed = i,
                                    save_dir1 = paste0(save_dir1, "/rand/"),
                                    file_index = i,
                                    randinds2 = randinds2)
      cat(i, " ")
    }
  }

  if(plot_real_vs_rand){

    # Load random indices:
    randpvals <- lapply(list.files(paste0(save_dir1,"/rand"), full.names = T), \(x) fread(x, data.table = F))

    ## Calculate # of mods sig in at least 1 subclass:

    # Calculate percentage of significant modules (real modules)
    realdf <- lapply(list(pvals), \(j){
      vals_bonf <- lapply(1:nrow(j),\(x) sum(j[x, ] < 0.05/(nrow(j)*ncol(j)))) |> unlist() 
      jq <- qval_redund(unlist(j))$qvalues |>
        matrix(ncol = ncol(j))
      vals_fdr <- lapply(1:nrow(j), \(x) sum(jq[x, ] < 0.05)) |> unlist()

      out <- data.frame("type" = "real",
                        "cutoff" = c("Nominal", "FDR", "Bonf"),
                        "pctg" = c(sum(apply(j, 1, \(x) sum(x < 0.05)) > 0) / 1023,
                                  sum(vals_fdr > 0) / 1023,
                                  sum(vals_bonf > 0) / 1023
                                  )
                        )
    }) |> unlist(recursive = F) |> as.data.frame()

    # Calculate percentage of significant modules (rand modules)
    randdf <- lapply(randpvals, \(j){
      vals_bonf <- lapply(1:nrow(j),\(x) sum(j[x, ] < 0.05/(nrow(j)*ncol(j)))) |> unlist() 
      jq <- qval_redund(unlist(j))$qvalues |>
        matrix(ncol = ncol(j))
      vals_fdr <- lapply(1:nrow(j), \(x) sum(jq[x, ] < 0.05)) |> unlist()
      out <- data.frame("type" = "rand",
                        "cutoff" = c("Nominal", "FDR", "Bonf"),
                        "pctg" = c(sum(apply(j, 1, \(x) sum(x < 0.05)) > 0) / 1023,
                                  sum(vals_fdr > 0) / 1023,
                                  sum(vals_bonf > 0) / 1023
                                  )
                        )
    }) |> do.call(what = "rbind")

    plotdf <- rbind(realdf, randdf)
    fwrite(plotdf, file = paste0(save_dir1, "/pctg_of_sig_modules_real_vs_rand.csv"))


    p <- plotdf |>
      mutate(pctg = pctg*100,
            cutoff = factor(cutoff, levels = unique(cutoff)),
            type = case_match(type, "real" ~ "Real network", "rand" ~ "Random networks\n(n = 1000)", .default = type),
            type = factor(type, levels = unique(type))) |>
      ggplot(aes(x = cutoff, y = pctg, fill = type, color = type)) +
        theme_classic() + 
        geom_bar(position = "dodge", stat = "summary", fun = "mean", alpha = 0.4, color = NA) +
        geom_point(position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.9), size = 0.2) +
        ylim(0, 100) +
        labs(x = "", y = "% of mods significant in >=1 subclass",
            title = "# of significant modules in at least 1 subclass using dCoPA",
            subtitle = "Real (n = 1023) vs random (n = 1023 x 1000) modules; 24 subclasses",
            caption = paste0("Significance test conducted using one-sided Wilcoxon\n", caption1)) +
        theme(legend.title = element_blank(),
              legend.position = "bottom") +
      ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", method.args = list(alternative = "less"), label.y = 80) + 
      ggpubr::stat_compare_means(method = "wilcox.test", aes(label = paste0("p = ", after_stat(p.format))),method.args = list(alternative = "less"), label.y = 90) +
      #ggpubr::stat_compare_means(method = "wilcox.test", aes(label = paste0("p = ", c(0.042, 0.012, 0.012))), label.y = 90) +
      guides(color = "none")

    #ggsave(p, file = paste0(save_dir1, "/pctg_of_sig_modules_real_vs_rand.png"), width = 6, height = 4)

    #################
    # ii) the % of all possible subclasses that show significant DE over all modules (i.e., number of significant subclasses / (# subclasses X # modules))
    #################
    # Calculate percentage of significant subclasses (real modules)
    realdf_ct <- lapply(list(pvals), \(j){
      vals_bonf <- lapply(1:nrow(j),\(x) sum(j[x, ] < 0.05/(nrow(j)*ncol(j)))) |> unlist() 
      jq <- qval_redund(unlist(j))$qvalues |>
        matrix(ncol = ncol(j))
      vals_fdr <- lapply(1:nrow(j), \(x) sum(jq[x, ] < 0.05)) |> unlist()
      out <- data.frame("type" = "real",
                        "cutoff" = c("Nominal", "FDR", "Bonf"),
                        "pctg" = c(sum(apply(j, 1, \(x) sum(x < 0.05))) / (1023 * 24),
                                  sum(vals_fdr) / (1023 * 24),
                                  sum(vals_bonf) / (1023 * 24)
                                  )
                        )
    }) |> unlist(recursive = F) |> as.data.frame()

    randdf_ct <- lapply(randpvals, \(j){
      vals_bonf <- lapply(1:nrow(j),\(x) sum(j[x, ] < 0.05/(nrow(j)*ncol(j)))) |> unlist() 
      jq <- qval_redund(unlist(j))$qvalues |>
        matrix(ncol = ncol(j))
      vals_fdr <- lapply(1:nrow(j), \(x) sum(jq[x, ] < 0.05)) |> unlist()
      out <- data.frame("type" = "rand",
                        "cutoff" = c("Nominal", "FDR", "Bonf"),
                        "pctg" = c(sum(apply(j, 1, \(x) sum(x < 0.05))) / (1023 * 24),                              
                                  sum(vals_fdr) / (1023 * 24),
                                  sum(vals_bonf) / (1023 * 24)
                                  )
                        )
    }) |> do.call(what = "rbind")

    plotdf_ct <- rbind(realdf_ct, randdf_ct)
    fwrite(plotdf_ct, file = paste0(save_dir1, "/pctg_of_sig_subclasses_real_vs_rand.csv"))


    p_ct <- plotdf_ct |>
      mutate(pctg = pctg*100,
            cutoff = factor(cutoff, levels = unique(cutoff)),
            type = case_match(type, "real" ~ "Real network", "rand" ~ "Random networks\n(n = 1000)", .default = type),
            type = factor(type, levels = unique(type))) |>
      ggplot(aes(x = cutoff, y = pctg, fill = type, color = type)) +
        theme_classic() + 
        geom_bar(position = "dodge", stat = "summary", fun = "mean", alpha = 0.4, color = NA) +
        geom_point(position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.9), size = 0.2) +
        ylim(0, 100) +
        labs(x = "", y = "% of subclasses that are significant",
            title = "# of significant subclasses over all modules using dCoPA",
            subtitle = "Real (n = 1023) vs random (n = 1023 x 1000) modules; 24 subclasses",
            caption = paste0("Significance test conducted using one-sided Wilcoxon\n", caption1)) +
        theme(legend.title = element_blank(),
              legend.position = "bottom") +
      ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif",,method.args = list(alternative = "less"), label.y = 80) + 
      ggpubr::stat_compare_means(method = "wilcox.test", aes(label = paste0("p = ", after_stat(p.format))),,method.args = list(alternative = "less"), label.y = 90) +
      guides(color = "none")

    #ggsave(p_ct, file = paste0(save_dir1, "/pctg_of_sig_subclasses_real_vs_rand.png"), width = 6, height = 4)

    pall <- cowplot::plot_grid(p, p_ct, nrow = 1)
    ggsave(pall, file = paste0(save_dir1,"/combined.png"), width = 12, height = 4)
  }

  #################
  # iii) the dCOPA ‘scorecard’ that reports the # of modules that are significant for each subclass (maybe as a dot plot). How long would it take to do this for n = 10, 100, or 1000 sets of random modules?
  #################

  ## Create dotplot
  # Load directionality 
  mod_direc <- fread(data.table = F, file = paste0(save_dir1, "/mod_direction.csv"))

  # Create vector of p-val cutoffs
  jq <- qval_redund(unlist(pvals))
  cutoffsp <- c(0.05, 
                jq$pvalues[which(jq$qvalues == max(jq$qvalues[jq$qvalues < 0.05]))] |> unique(),
                0.5/(24*1023))
  cutoffs <- -log10(cutoffsp)

  # Log transform pvalues (set 0 value to minimum excluding 0)
  lab1 <- c(paste0("Significant (FDR, higher in ", splits[1], " vs ", splits[2], ")"),
            paste0("Significant (FDR, lower in ", splits[1], " vs ", splits[2], ")"))
  pvals1 <- pvals
  pvals1[pvals1 == 0] <- min(pvals1[pvals1 != 0])
  pvals_log <- apply(pvals1, 2, \(x) -log10(x)) |> 
    as.data.frame() |>
    mutate(module = factor(1:1023),
          module = factor(module, levels = rev(module))) |>
    pivot_longer(!module, names_to = "subclass", values_to = "pvalue") |>
    mutate(dir = unlist(as.data.frame(t(mod_direc))),
          type = "real",
          pcut = case_when(
            #pvalue > cutoffs[1] & pvalue < cutoffs[2] ~ "Nominal",
            pvalue > cutoffs[2] & pvalue < cutoffs[3] & dir > 0 ~ lab1[1],
            pvalue > cutoffs[2] & pvalue < cutoffs[3] & dir < 0 ~ lab1[2],
          # pvalue > cutoffs[3] ~ "Bonf",
            .default = "Not significant"
    ))

  # pdot <- pvals_log |>
  #   ggplot(aes(x = subclass, y = module)) +
  #     theme_minimal() + 
  #     geom_point(aes(fill = pcut), color = "black", pch = 21) +
  #     theme(text = element_text(size = 6),
  #           axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
  #           axis.text.y = element_text(size = 6), 
  #           legend.position = "bottom",
  #           legend.title = element_text(hjust = 0.5),
  #           legend.spacing.y = unit(0, "mm")) +
  #     labs(x = "", y = "", fill = expression(paste("P-value (", -log[10], ")"))) +
  #     #scale_fill_manual(values = c("Nominal" = "green", "FDR" = "red", "Bonf" = "yellow", "Not significant" = "white")) +
  #     scale_fill_manual(values = c(lab1[2] = "green", lab1[1] = "red", "Not significant" = "white")) +
  #     guides(fill = guide_legend(title.position = "top", ncol = 1)) 

  # ggsave(pdot, file = paste0(save_dir1, "/dcopa_scorecard.pdf"), width = 3, height = 80, bg = "white", limitsize=F)

  # Separate plots per direction

  pvals_log_hi <- apply(pvals1, 2, \(x) -log10(x)) |> 
    as.data.frame() |>
    mutate(module = factor(1:1023),
          module = factor(module, levels = rev(module))) |>
    pivot_longer(!module, names_to = "subclass", values_to = "pvalue") |>
    mutate(dir = unlist(as.data.frame(t(mod_direc))),
          type = "real",
          pcut = case_when(
            #pvalue > cutoffs[1] & pvalue < cutoffs[2] ~ "Nominal",
            pvalue > cutoffs[2] & pvalue < cutoffs[3] & dir > 0 ~ lab1[1],
          # pvalue > cutoffs[3] ~ "Bonf",
            .default = "Not significant"
    ))

  scale_val <- c("red", "white")
  names(scale_val) <- c(lab1[1], "Not significant")
  pdothi <- pvals_log_hi |>
    ggplot(aes(x = subclass, y = module)) +
      theme_minimal() + 
      geom_point(aes(fill = pcut), color = "black", pch = 21) +
      theme(text = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 6), 
            legend.position = "bottom",
            legend.title = element_text(hjust = 0.5),
            legend.spacing.y = unit(0, "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm"))) +
      labs(x = "", y = "", fill = expression(paste("P-value (", -log[10], ")")), 
          title = paste0("Higher in ", splits[1], " vs ", splits[2]),
          caption = caption1) +
      scale_fill_manual(values = scale_val) +
      guides(fill = guide_legend(title.position = "top", nrow = 1)) 

  #ggsave(pdothi, file = paste0(save_dir1, "/dcopa_scorecard_higherInAD.pdf"), width = 3, height = 80, bg = "white", limitsize=F)
 

  pvals_log_lo <- apply(pvals1, 2, \(x) -log10(x)) |> 
    as.data.frame() |>
    mutate(module = factor(1:1023),
          module = factor(module, levels = rev(module))) |>
    pivot_longer(!module, names_to = "subclass", values_to = "pvalue") |>
    mutate(dir = unlist(as.data.frame(t(mod_direc))),
          type = "real",
          pcut = case_when(
            #pvalue > cutoffs[1] & pvalue < cutoffs[2] ~ "Nominal",
            pvalue > cutoffs[2] & pvalue < cutoffs[3] & dir < 0 ~ lab1[2],
          # pvalue > cutoffs[3] ~ "Bonf",
            .default = "Not significant"
    ))

  scale_val_lo <- c("green", "white")
  names(scale_val_lo) <- c(lab1[2], "Not significant")
  pdotlo <- pvals_log_lo |>
    ggplot(aes(x = subclass, y = module)) +
      theme_minimal() + 
      geom_point(aes(fill = pcut), color = "black", pch = 21) +
      theme(text = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 6), 
            legend.position = "bottom",
            legend.title = element_text(hjust = 0.5),
            legend.spacing.y = unit(0, "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm"))) +
      labs(x = "", y = "", fill = expression(paste("P-value (", -log[10], ")")), 
          title = paste0("Lower in ", splits[1], " vs ", splits[2]),
          caption = caption1) +
      scale_fill_manual(values = scale_val_lo) +
      guides(fill = guide_legend(title.position = "top", nrow = 1)) 

  pdotall <- cowplot::plot_grid(pdothi, pdotlo, ncol = 2)
  ggsave(pdotall, file = paste0(save_dir1, "/dcopa_scorecard.pdf"), width = 6, height = 80, bg = "white", limitsize=F)


  ## Create summary dotplots (# of significant mods per subclass)
  # Add colors for cell class
  sn_anno_subclass <- fread(data.table = F, file = file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv")) |> 
    select(RNA.Class, RNA.Subclass) |>
    filter(!duplicated(RNA.Subclass), 
           RNA.Subclass %in% colnames(pvals1)) |>
    mutate(Class = case_match(RNA.Class, "Exc" ~ "Glutamatergic", "Inh" ~ "GABAergic", .default = "Non-neuronal")) |>
    rename(Subclass = RNA.Subclass) |>
    arrange(factor(Class, levels = c("Glutamatergic", "GABAergic", "Non-neuronal")))

  hi_sum <- pvals_log_hi |>
    filter(pcut != "Not significant") |>
    group_by(pcut, subclass) |>
    summarise(num_sig = n()) |>
    left_join(sn_anno_subclass, by = join_by(subclass == Subclass)) |>
    mutate(type = paste0("Higher in ", splits[1], " vs ", splits[2])) |>
    as.data.frame()

  cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

  # lower in AD
  lo_sum <- pvals_log_lo |>
    filter(pcut != "Not significant") |>
    group_by(pcut, subclass) |>
    summarise(num_sig = n()) |>
    left_join(sn_anno_subclass, by = join_by(subclass == Subclass)) |>
    mutate(type = paste0("Lower in ", splits[1], " vs ", splits[2])) |>
    as.data.frame()

  psumall <- rbind(hi_sum, lo_sum) |>
    arrange(factor(subclass, levels = sn_anno_subclass$Subclass)) |>
    mutate(pcut = "",
          type = factor(type, levels = rev(unique(type))),
          subclass = factor(subclass, levels = unique(subclass))) |>
    ggplot(aes(x = subclass, y = type, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.2, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 6), 
            strip.text = element_text(size = 7),
            legend.position = "bottom",
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(x = "", y = "",
          title = "# of significant dCoPA mods",
          subtitle = "1023 total modules, FDR cutoff",
          caption = caption1) +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(scales::pretty_breaks(n = 3)) 

  ggsave(psumall, file = paste0(save_dir1, "/dcopa_scorecard_summary.pdf"), width = 5, height = 2.2, bg = "white", limitsize=F)
  ggsave(psumall, file = paste0(save_dir1, "/dcopa_scorecard_summary.png"), width = 5, height = 2.2, bg = "white", limitsize=F)

  ## Additional exploratory plots
  # Size of the module (# genes) vs # of significant subclasses identified by dCoPA?  
  # What about PC1 VE and the # of significant subclasses?  
  # Could look at all significant subclasses or stratify by those that are up in AD / down in AD. 

  # # Calc % VE
  # bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
  mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
  mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
    (\(x) tapply(x[, 2], x[, 3], list))()
  mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
  these_mods <- which(mod_bc_lengths > 3)
  mod_seed_lengths <- lapply(mod_seed,length) |> unlist()
  # expr_t <- t(bulk_expr[,3:ncol(bulk_expr)])
  # colnames(expr_t) <- bulk_expr[,2]
  # expr_z <- apply(expr_t, 2, function(x) (x-mean(x))/sd(x))
  # r2_vec_bc_list <- list()
  # for(j in seq_along(mod_seed)){
  #   expr_plot <- as.data.frame(expr_z[ ,colnames(expr_z) %in% mod_seed[[j]]])
  #   pc1_bc <-  prcomp(expr_plot[apply(expr_plot,2,var)>0], scale=T)$x[,1]
  #   r2_vec_bc <- apply(expr_plot[apply(expr_plot,2,var)>0],2,function(x) summary(lm(x~pc1_bc))$r.squared)
  #   r2_vec_bc_list[[j]] <- r2_vec_bc
  #   cat(j, " ")
  # }
  # vevec <- lapply(r2_vec_bc_list, mean) |> unlist()
  # qsave(vevec, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/PC1VE_1023.qs"))
  vevec <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/PC1VE_1023.qs"))


  desc_mat <- data.frame(
    "mod_seed_size" = mod_seed_lengths,
    "mod_bc_size" = mod_bc_lengths,
    "VE" = vevec
  )

  # Filter to 1023
  desc_mat <- desc_mat[these_mods, ]

  # Add significance info
  desc_mat$sig_count_all <- apply(pvals1, 1, \(x) sum(x < cutoffsp[2]))
  desc_mat$sig_count_hi <- apply(apply(pvals1, 2, \(x) x < cutoffsp[2]) & apply(mod_direc, 2, \(x) x > 0), 1, sum)
  desc_mat$sig_count_lo <- apply(apply(pvals1, 2, \(x) x < cutoffsp[2]) & apply(mod_direc, 2, \(x) x < 0), 1, sum)

  desc_cors <- cor(desc_mat) |>
    apply(2, \(x) sprintf("%.3g",  x))

  # cor(desc_mat)
  #               mod_seed_size mod_bc_size         VE sig_count_all sig_count_hi
  # mod_seed_size  1.0000000000  0.15588867 0.13264213     0.1430074 0.0009364814
  # mod_bc_size    0.1558886713  1.00000000 0.16427111     0.1853054 0.0308388753
  # VE             0.1326421283  0.16427111 1.00000000     0.3601409 0.0268212672
  # sig_count_all  0.1430073631  0.18530537 0.36014085     1.0000000 0.2096452486
  # sig_count_hi   0.0009364814  0.03083888 0.02682127     0.2096452 1.0000000000
  # sig_count_lo   0.1460587793  0.18320732 0.36282797     0.9798480 0.0101149596
  #               sig_count_lo
  # mod_seed_size   0.14605878
  # mod_bc_size     0.18320732
  # VE              0.36282797
  # sig_count_all   0.97984805
  # sig_count_hi    0.01011496
  # sig_count_lo    1.00000000

  dplotlist <- list()

  dplotlist[[1]] <- ggplot(desc_mat, aes(x = sig_count_all, y = mod_seed_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = "# of significant subclasses (all)", y = "Module size (seed)", title = paste0("r = ", desc_cors[1, 4])) +
    geom_smooth() 
  dplotlist[[2]] <- ggplot(desc_mat, aes(x = sig_count_all, y = mod_bc_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = "# of significant subclasses (all)", y = "Module size\n(topmodposbc)", title = paste0("r = ", desc_cors[2, 4])) +
    geom_smooth()
  dplotlist[[3]] <- ggplot(desc_mat, aes(x = sig_count_all, y = mod_bc_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = "# of significant subclasses (all)", y = "PC1 %VE", title = paste0("r = ", desc_cors[3, 4])) +
    geom_smooth()
  dplotlist[[4]] <- ggplot(desc_mat, aes(x = sig_count_hi, y = mod_seed_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = paste0("# of significant subclasses\n(higher in ",splits[1], " vs ", splits[2], ")"), y = "Module size (seed)", title = paste0("r = ", desc_cors[1, 5])) +
    geom_smooth()
  dplotlist[[5]] <- ggplot(desc_mat, aes(x = sig_count_hi, y = mod_bc_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = paste0("# of significant subclasses\n(higher in ",splits[1], " vs ", splits[2], ")"), y = "Module size\n(topmodposbc)", title = paste0("r = ", desc_cors[2, 5])) +
    geom_smooth()
  dplotlist[[6]] <- ggplot(desc_mat, aes(x = sig_count_hi, y = mod_bc_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = paste0("# of significant subclasses\n(higher in ",splits[1], " vs ", splits[2], ")"), y = "PC1 %VE", title = paste0("r = ", desc_cors[3, 5])) +
    geom_smooth()
  dplotlist[[7]] <- ggplot(desc_mat, aes(x = sig_count_lo, y = mod_seed_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = paste0("# of significant subclasses\n(lower in ",splits[1], " vs ", splits[2], ")"), y = "Module size (seed)", title = paste0("r = ", desc_cors[1, 6])) +
    geom_smooth()
  dplotlist[[8]] <- ggplot(desc_mat, aes(x = sig_count_lo, y = mod_bc_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = paste0("# of significant subclasses\n(lower in ",splits[1], " vs ", splits[2], ")"), y = "Module size\n(topmodposbc)",title = paste0("r = ", desc_cors[2, 6])) +
    geom_smooth()
  dplotlist[[9]] <- ggplot(desc_mat, aes(x = sig_count_lo, y = mod_bc_size)) +
    theme_classic() + 
    geom_point() +
    labs(x = paste0("# of significant subclasses\n(lower in ",splits[1], " vs ", splits[2], ")"), y = "PC1 %VE", title = paste0("r = ", desc_cors[3, 6])) +
    geom_smooth()

  dplotall <- cowplot::plot_grid(plotlist = dplotlist, nrow = 3, ncol = 3)
  ggsave(dplotall, file = paste0(save_dir1, "/summary_plots.png"), width = 11, height = 7)

}

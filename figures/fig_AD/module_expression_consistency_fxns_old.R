library(ggpubr)
library(showtext)
showtext_auto()

### Plotting function
plotmods <- function(subdiffallabs, 
                     subdiffallper, 
                     adconmeans,
                     indsmat,
                     save_path,
                     file_name,
                     all_plots,
                     sea_means
                     ){

  inds <- unique(indsmat$mod)

  cols <- RColorBrewer::brewer.pal(6, "Paired")
  cols2 <- RColorBrewer::brewer.pal(10, "Spectral")

  out_plot <- lapply(seq_along(inds), \(i){
    j <- which(these_mods == inds[i])
    newtheme <- theme_light() + theme(axis.text.x = element_text(size = 4), 
                                    axis.text.y = element_text(size = 4), 
                                    axis.title.x = element_text(size = 4), 
                                    axis.title.y = element_text(size = 4), 
                                    legend.text = element_text(size = 4),
                                    plot.subtitle = element_blank(),
                                    legend.key.size = unit(0.2, "cm"), 
                                    legend.title = element_blank())
    plots <- lapply(all_plots, \(x) x[[j]] + newtheme)

    plots[[1]]$layers[[1]] <- NULL
    plots[[1]] <- plots[[1]] + 
        geom_line(linewidth = 0.2) +
        labs(x = "Sample") + 
        theme(axis.title.x = element_text(size = 4), 
            axis.title.y = element_text(size = 4),
            legend.text = element_text(size = 4, margin = margin(0, 0, 0, -0.1)), 
            legend.box.margin = margin(0, 0, 0, -7),
            legend.position = "right",
            plot.title = element_text(size = 6, hjust = 0.5, margin = margin(0, 0, -0.1, 0))) +
        scale_color_manual(values = cols2) +
        ggtitle(paste0("Module ", inds[i])) +
        guides(color=guide_legend(ncol=2,theme = theme(legend.byrow = F)))
        

    plots[[2]] <- all_plots[[2]][[j]]$data |> 
        dplyr::filter(pval>0) |>
        dplyr::mutate(SetName=factor(SetName, levels=rev(unique(SetName)))) |>
        ggplot(aes(x = SetName, y = pval)) +
            theme_light() +
            geom_bar(stat="identity") +
            theme(axis.text.x = element_text(size=6,hjust=1,vjust = 0.5),
                #axis.title.x = element_text(size=6),
                plot.title = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size=3),
                axis.text.y = element_text(size=3),
                plot.margin = margin(0,0,0,-13)) +
            labs(x="", y=bquote(-log[10](p-val))) #+
            #scale_x_discrete(limits = rev(levels(SetName))) +
            #geom_hline(yintercept = gsea_cutFDR, color = "red")
    plots[[2]]$layers[[2]] <- all_plots[[2]][[j]]$layers[[2]]
    plots[[2]] <- plots[[2]] + coord_flip()
    #p <- plot_grid(plotlist = plots, ncol = 4, nrow = 1, align = "h", axis = "bt", rel_widths=c(0.5, 0.85, 0.75, 1, 0.8))
    #return(p)

    plots[[3]] <- adconmeans |> 
    dplyr::filter(mod == inds[i]) |>
    ggplot(aes(x = subclass, y = mean, fill = type)) + 
        theme_classic() + 
        geom_boxplot(outlier.size = 0.1, linewidth = 0.2) +
        theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=4),
              axis.text.y = element_text(size=4),
              legend.position = "right",
              legend.margin = margin(0, 0, 0, -10),
              legend.text = element_text(size = 4, margin = margin(0, 0, 0, -0.1)),
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              plot.title = element_blank()) +
        labs(y = "Mean expression")


    #highlight_ct <- indsmat$subclass[i]
    
    plots[[4]] <- subdiffallabs |> 
      dplyr::filter(modno == inds[i]) |>
      mutate(highlight = factor(ifelse(highlight, "Significant", "Not significant"), levels = c("Significant", "Not significant"))) |>
      #mutate(highlight = subclass == highlight_ct) |>
      ggplot(aes(x = subclass, y = meanpd, fill = highlight)) +
        theme_classic() + 
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=4),
            axis.text.y = element_text(size=4),
           # legend.margin = margin(0, 0, 0, -10),
         #   legend.text = element_text(size = 4, margin = margin(0, 0, 0, -0.1)),
          #  legend.key.size = unit(0.2, "cm"),
            legend.position = "none",
          #  legend.title = element_blank(),
            plot.title = element_blank()) + 
        labs(y = paste0("Mean change in expr\n(", names(sea_means)[1], " minus ", names(sea_means)[2], ")"))
    
    plots[[5]] <-  subdiffallper |> 
        dplyr::filter(modno == inds[i]) |>
        mutate(highlight = factor(ifelse(highlight, "Significant", "Not significant"), levels = c("Significant", "Not significant"))) |>
        #mutate(highlight = subclass == highlight_ct) |>
        ggplot(aes(x = subclass, y = meanpd, fill = highlight)) +
            theme_classic() + 
            geom_bar(stat = "identity") +
            theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size=4),
                axis.text.y = element_text(size=4),
                legend.margin = margin(0, 0, 0, -10),
                legend.text = element_text(size = 4, margin = margin(0, 0, 0, -0.1)),
                legend.key.size = unit(0.2, "cm"),
                legend.position = "right",
                legend.title = element_blank(),
                plot.title = element_blank()) + 
            labs(y = paste0("Mean % change\nin expr\n(", names(sea_means)[1], " minus ", names(sea_means)[2], ")"))

    plots[[6]] <- mapply(\(x, y){
      highlight_cts <- subdiffallabs |>
        dplyr::filter(highlight,
                      modno == inds[i]) |>
        dplyr::pull(subclass)
      df <- x[rownames(x) %in% mod_bc[[inds[i]]], colnames(x) %in% highlight_cts, drop = F] 
      out <- data.frame("type" = y, df) |>
        pivot_longer(!type, names_to = "ct", values_to = "vals")
      return(out)
    }, sea_means[1:2], names(sea_means)[1:2], SIMPLIFY = F) |>
      do.call(what = "rbind") |>
      ggpaired(x = "type", 
               y = "vals", 
               color = "type", 
               line.color = "gray", 
               line.size = 0.4,
               point.size = 0.4,
               notch = T) +
        facet_wrap(~ct, nrow = 1) + 
        labs(x = "", y = "Mean expression") +
        theme(legend.position = "none", 
            axis.text.x = element_text(size = 4),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 4),
            axis.text.y = element_text(size = 4),
            strip.text = element_text(size = 4, margin = margin(-5, 0, -5, 0)))  
    return(plots)
  })

  # Plot individual snapshots
  save_indiv <- file.path(save_path, "indiv")
  if(!dir.exists(save_indiv)){dir.create(save_indiv, recursive = T)}
  outalllist <- mapply(\(x, y){
    # p <- cowplot::plot_grid(plotlist = x, ncol = 6, align = "h", axis = "b", rel_widths = c(0.6, 0.7, 0.8, 0.8, 0.8, 0.4)) +
    #   coord_cartesian(clip = "off")
    p1 <- cowplot::plot_grid(plotlist = x[1:3], ncol = 3, align = "h", axis = "b", rel_widths = c(0.6, 0.7, 0.8)) +
      coord_cartesian(clip = "off")
    p2 <- cowplot::plot_grid(plotlist = x[4:6], ncol = 3, align = "h", axis = "b", rel_widths = c(0.8, 1, 1)) +
      coord_cartesian(clip = "off")
    p <- cowplot::plot_grid(p1, p2, nrow = 2) 
    #p <- cowplot::plot_grid(plotlist = x, ncol = 3, nrow = 2, align = "h", axis = "b")
    ggsave(p, file = file.path(save_indiv, paste0(sprintf("%04d", as.numeric(y)), ".pdf")),  height = 2, width = 7, bg = "white", limitsize = F)
    ggsave(p, file = file.path(save_indiv, paste0(sprintf("%04d", as.numeric(y)), ".svg")),  height = 2, width = 7, bg = "white", limitsize = F)
   # ggsave(p, file = file.path(save_indiv, paste0(y, ".png")),  height = 2, width = 9, bg = "white", limitsize = F)
  }, out_plot, inds, SIMPLIFY = F)

  # out_plot2 <- lapply(c(1,2,3,4,5,6), \(x){
  #   plist <- lapply(out_plot, \(y) y[[x]])
  #   if(x != 2){
  #       return(plot_grid(plotlist = plist, nrow = length(plist), align = "v", axis = "rl"))
  #   } else {
  #       return(plot_grid(plotlist = plist, nrow = length(plist)))
  #   }
  #   cat(x, "\n")
  # })

  # # Plot all shapshots together
  # outall2 <- cowplot::plot_grid(plotlist = out_plot2, ncol = 6, align = "h", axis = "bt", rel_widths = c(0.6, 0.7, 0.8, 0.8, 0.8, 0.4))
  # ggsave(outall2, file = file.path(save_path, paste0(file_name, ".png")),  height = length(inds) * 1, width = 15, bg = "white", limitsize = F)
  # ggsave(outall2, file = file.path(save_path, paste0(file_name, ".pdf")),  height = length(inds) * 1, width = 15, bg = "white", limitsize = F)

}
 
# Wrapper
find_cons_mods <- function(save_path,
                           copa_dir1,
                           sn_summary_path = NULL,
                           return_obj = F,
                           subset = NULL, # vector of subset of modules to plot
                           plot = T,
                           comp_names = NULL
                           ){
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = T)
  }

  # Load objects
  euc_dist <- list("all"=qread(paste0(copa_dir1, "/euclidean_distances/euclidean_sigmods_all.qs")),
                   "pos"=qread(paste0(copa_dir1, "/euclidean_distances/euclidean_sigmods_positive.qs")),
                   "neg"=qread(paste0(copa_dir1, "/euclidean_distances/euclidean_sigmods_negative.qs")))
  if(is.null(sn_summary_path)){
    sea_means <- qread(paste0(copa_dir1, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
  } else {
    sea_means <- qread(paste0(sn_summary_path, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
  }
  if(!is.null(comp_names)){
    names(sea_means) <- comp_names
  }
  sea_means[[2]] <- sea_means[[2]][match(rownames(sea_means[[1]]), rownames(sea_means[[2]])), ]
  sea_means <- lapply(sea_means, as.data.frame)
  sub_diff <- sea_means[[1]] - sea_means[[2]] |> as.data.frame()
  sub_diff_per <- ((sea_means[[1]] - sea_means[[2]]) / sea_means[[1]]) 
  seadist <- fread(data.table = F, file = paste0(copa_dir1, "/euclidean_distances/euclidean_distances_all_Subclass.csv"))
  rownames(seadist) <- these_mods
  all_plots <- list(qread(file.path(copa_dir1, "sn_proj_objects", "expr_object_bc.qs")), 
                    qread(file.path(copa_dir1, "sn_proj_objects", "gsea_object.qs")))

  # Find consistent modules (pos)
  copacohpos_complete <- mapply(\(subclass, subclass_name){
      outvec <- lapply(subclass, \(mod){
          temp <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], colnames(sub_diff) == subclass_name]
          return(sum(temp > 0) / length(temp))
      }) |> unlist() 
      if(length(outvec) > 0){
        return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
      } else {
        return(data.frame("subclass" = subclass_name, "pcnt" = NA))
      }
  }, euc_dist$pos[[3]], names(euc_dist$pos[[3]]), SIMPLIFY = F) |>
    lapply(\(x){
      x[x$pcnt == 1, ] |>
      rownames_to_column(var = "mod")
    }) |> do.call(what = "rbind") |>
      dplyr::filter(!is.na(subclass))

  # Find consistent modules (neg)
  copacohneg_complete <- mapply(\(subclass, subclass_name){
    outvec <- lapply(subclass, \(mod){
      temp <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], colnames(sub_diff) == subclass_name]
      return(sum(temp < 0) / length(temp))
    }) |> unlist() 
    if(length(outvec) > 0){
      return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
    } else {
      return(data.frame("subclass" = subclass_name, "pcnt" = NA))
    }
  }, euc_dist$neg[[3]], names(euc_dist$neg[[3]]), SIMPLIFY = F) |>
    lapply(\(x){
      x[x$pcnt == 1, ] |>
        rownames_to_column(var = "mod")
    }) |> do.call(what = "rbind") |>
      dplyr::filter(!is.na(subclass))

  if(nrow(copacohneg_complete) == 0 | nrow(copacohpos_complete) == 0){
    stop("No modules remain after filtering")
  }

  # pos_copa_unique <- lapply(1:nrow(copacohpos_complete), \(i){
  #   temp <- names(euc_dist$pos[[3]])[unlist(lapply(euc_dist$pos[[3]], \(x) copacohpos_complete$mod[i] %in% x))]
  #   return(temp)
  # })
  
  copacohpos_complete <- copacohpos_complete |>
    mutate(max_AD = lapply(mod, \(modind){ # Add info about subclass with highest expression in AD
      mean1 <- sea_means[[1]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
      return(names(mean1)[which.max(mean1)])
    }) |> unlist()) |>
    mutate(max_con = lapply(mod, \(modind){ # Add info about subclass with highest expression in con
      mean1 <- sea_means[[2]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
      return(names(mean1)[which.max(mean1)])
    }) |> unlist()) |>
    mutate(dist_all = lapply(mod, \(modind){ # Add info about total euclidean distance (over all subclasses) 
      seadist$all[which(these_mods %in% modind)]
    }) |> unlist()) |>
    mutate(dist_per = lapply(mod, \(modind){ # Add info about subclass with highest percentage diff
      tempdf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[modind]], ] |> colMeans() |> unlist()
      return(names(tempdf)[which.max(tempdf)])
    }) |> unlist()) |>
    arrange(desc(dist_all)) |> 
    arrange(as.numeric(mod))

  # filt_these <- mapply(\(x, y){
  #   x %in% y
  # }, copacohpos_complete$subclass, pos_copa_unique, SIMPLIFY = F) |> unlist()

  # copacohpos_complete <- copacohpos_complete[filt_these, ] |>
  #   #dplyr::filter(subclass == max_AD | subclass == max_con) |>
  #   arrange(desc(dist_all)) |> 
  #   arrange(as.numeric(mod))
  
  
  # Repeat the above steps for negative (higher in con) modules
  
  # neg_copa_unique <- lapply(1:nrow(copacohneg_complete), \(i){
  #   temp <- names(euc_dist$neg[[3]])[unlist(lapply(euc_dist$neg[[3]], \(x) copacohneg_complete$mod[i] %in% x))]
  #   return(temp)
  # })

  copacohneg_complete <- copacohneg_complete |>
    mutate(
      max_AD = lapply(mod, \(modind){
        mean1 <- sea_means[[1]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> 
          apply(2, mean) |> 
          unlist()
        return(names(mean1)[which.max(mean1)])
      }) |> unlist(),
      max_con = lapply(mod, \(modind){
        mean1 <- sea_means[[2]][rownames(sea_means[[1]]) %in% mod_bc[[modind]], ] |> apply(2, mean) |> unlist()
        return(names(mean1)[which.max(mean1)])
      }) |> unlist(),
      dist_all = lapply(mod, \(modind){
        seadist$all[which(these_mods %in% modind)]
      }) |> unlist(),
      dist_per = lapply(mod, \(modind){ # Add info about subclass with highest percentage diff
        tempdf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[modind]], ] |> colMeans() |> unlist()
        return(names(tempdf)[which.min(tempdf)])
      }) |> unlist()
    ) |>
    # dplyr::filter((copa_unique == max_AD | copa_unique == max_con) & subclass == copa_unique) |>
    # arrange(desc(dist_all))
    arrange(desc(dist_all)) |> 
    arrange(as.numeric(mod))

  # filt_these <- mapply(\(x, y){
  #   x %in% y
  # }, copacohneg_complete$subclass, neg_copa_unique, SIMPLIFY = F) |> unlist()

  # copacohneg_complete <- copacohneg_complete[filt_these, ] |>
  #   #dplyr::filter(subclass == max_AD | subclass == max_con) |>
  #   arrange(desc(dist_all)) |> 
  #   arrange(as.numeric(mod))

  # Filter for input subset of mods if available
  if(!is.null(subset)){
    copacohpos_complete <- copacohpos_complete |>
      dplyr::filter(mod %in% subset)
    copacohneg_complete <- copacohneg_complete |>
      dplyr::filter(mod %in% subset)
  }

  # Gather percentage differences (pos)
  thesead2 <- unique(copacohpos_complete[,1]) |> as.numeric()
  sub_diff_mod <- lapply(thesead2, \(mod){
    moddf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[mod]], ] |>
      pivot_longer(everything(), names_to = "subclass", values_to = "pcnt_diff") |>
      mutate("modno" = mod) 
      return(moddf)
    }) |> do.call(what = "rbind") |>
      group_by(subclass, modno) |>
      summarise("meanpd" = mean(pcnt_diff)) |>
      mutate("highlight" = F)
  for(i in 1:nrow(sub_diff_mod)){
    if(sub_diff_mod$subclass[i] %in% copacohpos_complete$subclass[copacohpos_complete$mod == sub_diff_mod$modno[i]]){
      sub_diff_mod$highlight[i] <- T
    }
  }  

  # Gather absolute differences (pos)
  sub_diff_mod1 <- lapply(thesead2, \(mod){
    moddf <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], ] |>
      pivot_longer(everything(), names_to = "subclass", values_to = "diff") |>
      mutate("modno" = mod) 
    return(moddf)
  }) |> 
    do.call(what = "rbind") |>
    group_by(subclass, modno) |>
    summarise("meanpd" = mean(diff)) |>
    mutate("highlight" = F)
  for(i in 1:nrow(sub_diff_mod1)){
    if(sub_diff_mod1$subclass[i] %in% copacohpos_complete$subclass[copacohpos_complete$mod == sub_diff_mod1$modno[i]]){
      sub_diff_mod1$highlight[i] <- T
    }
  }  

  # Gather percentage differences (neg)
  theseadneg2 <- unique(copacohneg_complete[,1]) |> as.numeric()
  sub_diff_modneg <- lapply(theseadneg2, \(mod){
    moddf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[mod]], ] |>
      pivot_longer(everything(), names_to = "subclass", values_to = "pcnt_diff") |>
      mutate("modno" = mod) 
    return(moddf)
  }) |> do.call(what = "rbind")
  if(!is.null(sub_diff_modneg)){
    sub_diff_modneg <- sub_diff_modneg |>
      group_by(subclass, modno) |>
      summarise("meanpd" = mean(pcnt_diff)) |>
      mutate("highlight" = F)
  for(i in 1:nrow(sub_diff_modneg)){
      if(sub_diff_modneg$subclass[i] %in% copacohneg_complete$subclass[copacohneg_complete$mod == sub_diff_modneg$modno[i]]){
        sub_diff_modneg$highlight[i] <- T
      }
    }  
  } 

  # Gather absolute differences (neg)
  sub_diff_modneg1 <- lapply(theseadneg2, \(mod){
    moddf <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], ] |>
      pivot_longer(everything(), names_to = "subclass", values_to = "diff") |>
      mutate("modno" = mod) 
    return(moddf)
  }) |> do.call(what = "rbind") 
  if(!is.null(sub_diff_modneg)){
    sub_diff_modneg1 <- sub_diff_modneg1 |>
      group_by(subclass, modno) |>
      summarise("meanpd" = mean(diff)) |>
      mutate("highlight" = F)
    for(i in 1:nrow(sub_diff_modneg1)){
      if(sub_diff_modneg1$subclass[i] %in% copacohneg_complete$subclass[copacohneg_complete$mod == sub_diff_modneg1$modno[i]]){
        sub_diff_modneg1$highlight[i] <- T
      }
    }  
  }
  
  # Create dataframe for projection snapshot (boxplots)
  subdiffallper <- rbind(sub_diff_mod, sub_diff_modneg) |>
    mutate(meanpd = meanpd * 100)
  # subdiffallabs <- rbind(sub_diff_mod1, sub_diff_modneg1)
  adconmeans <- lapply(unique(subdiffallper[[2]]), \(mod){
    mapply(\(sea, name){
      sea[rownames(sea) %in% mod_bc[[mod]], ] |>
        pivot_longer(everything(), names_to = "subclass", values_to = "mean") |>
        mutate("type" = name)
    }, sea_means[1:2], names(sea_means)[1:2], SIMPLIFY = F) |>
      do.call(what = "rbind") |>
      mutate(mod = mod)
  }) |> do.call(what = "rbind")

  # Run plotting function
  if(plot){
    if(length(copacohpos_complete$mod) > 0){
      plotmods(subdiffallabs = sub_diff_mod1,
               subdiffallper = sub_diff_mod |> mutate(meanpd = meanpd * 100), 
               adconmeans = adconmeans,
               indsmat = copacohpos_complete,
               save_path = file.path(save_path, "pos"),
               file_name = paste0("modList_higherIn", names(sea_means)[1], "_filtered"),
               all_plots = all_plots,
               sea_means = sea_means
               )
    }
      
    if(length(copacohneg_complete$mod) > 0){
      plotmods(subdiffallabs = sub_diff_modneg1,
                subdiffallper = sub_diff_modneg |> mutate(meanpd = meanpd * 100), 
                adconmeans = adconmeans,
                indsmat = copacohneg_complete,
                save_path = file.path(save_path, "neg"),
                file_name = paste0("modList_higherIn", names(sea_means)[2], "_filtered"),
                all_plots = all_plots,
                sea_means = sea_means
                )
    }
  }

  if(return_obj){
    objs <- list(copacohpos_complete, copacohneg_complete)
    qsave(objs, file = file.path(save_path, "modlist_filtered.qs"))
    return(objs)
  }
}

# Plot differences between mean expression of all genes per subclass (barplots)
plot_allgenes_diff <- function(copa_dir1,
                               ytitle,
                               plottitle,
                               #cohort_names,
                               save_dir#,
                              # save_dir2
                               ){
    sea_means <- qread(paste0(copa_dir1, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
    sea_means[[2]] <- sea_means[[2]][match(rownames(sea_means[[1]]), rownames(sea_means[[2]])), ]
    sea_means <- lapply(sea_means, as.data.frame)
    sub_diff <- sea_means[[1]] - sea_means[[2]] |> as.data.frame()
    sub_diff_means <- apply(sub_diff, 2, mean)

    meantestsbonf <- mapply(\(ad, con){
        p <- wilcox.test(ad, con)$p.value
        return(p * 24)
    }, sea_means[[1]] |> as.list(), sea_means[[2]] |> as.list(), SIMPLIFY = F) |> unlist()

    plotdf <- data.frame("ct" = names(sub_diff_means), "meandiff" = sub_diff_means, "pval" = -log10(meantestsbonf))
    plotdf$star <- ""
    plotdf$star[plotdf$pval >= -log10(.05)] <- "*"
    plotdf$star[plotdf$pval >= -log10(.01)] <- "**"
    plotdf$star[plotdf$pval >= -log10(.001)] <- "***"
    plotdf$ypos = ifelse(plotdf$meandiff > 0, plotdf$meandiff, 0)

    p <- ggplot(plotdf, aes(x = ct, y = meandiff, fill = pval)) + 
    theme_classic() +
    geom_bar(stat = "identity") +
    geom_text(aes(label = star, y = ypos), vjust = -0.5) +
    theme(text = element_text(size = 18), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title.y = element_text(size = 16),
            plot.margin = margin(1, 1, 1, 1, "cm"),
            legend.title = element_text(hjust = 0.5, size = 12),
            plot.title = element_text(hjust = 0.5)) +
    labs(x = "", y = ytitle,
        fill = bquote(-log[10]~"p-val"),
        title = plottitle)
    ggsave(p, file = save_dir, width = 11, height = 6.5)

    # plotdflist <- list()
    # for(i in seq_along(sea_means[[1]])){
    #   plotdf <- data.frame("ct" = colnames(sea_means[[1]])[i],
    #                       "means" = c(sea_means[[1]][, i], sea_means[[2]][, i]),
    #                       "type" = c(rep(cohort_names[1], nrow(sea_means[[1]])), rep(cohort_names[2], nrow(sea_means[[1]]))),
    #                       "pval" = -log10(meantestsbonf[i])
    #                       )
    #   plotdf$star <- ""
    #   plotdf$star[plotdf$pval >= -log10(.05)] <- "*"
    #   plotdf$star[plotdf$pval >= -log10(.01)] <- "**"
    #   plotdf$star[plotdf$pval >= -log10(.001)] <- "***"
    #   plotdflist[[i]] <- plotdf
    # }
    # plotdfall <- do.call(rbind, plotdflist)
    # p <- ggplot(plotdfall, aes(x = ct, y = means, fill = type)) + 
    #     theme_classic() +
    #     geom_boxplot(outlier.size = 0) +
    #     #geom_text(aes(label = star), vjust = -0.5) +
    #     theme(text = element_text(size = 18), 
    #           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    #           axis.title.y = element_text(size = 16),
    #           plot.margin = margin(1, 1, 1, 1, "cm"),
    #           legend.title = element_blank(),
    #           plot.title = element_text(hjust = 0.5)) +
    #     labs(x = "", y = ytitle,
    #         fill = bquote(-log[10]~"p-val"),
    #         title = plottitle)
    # ggsave(p, file = save_dir2, width = 11, height = 6.5)
    
}

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
                     sea_means,
                     header_name
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
        geom_boxplot(outlier.size = 0.05, linewidth = 0.2) +
        theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=4),
              axis.text.y = element_text(size=4),
              legend.position = "right",
              legend.margin = margin(0, 0, 0, -10),
              legend.text = element_text(size = 4, margin = margin(0, 0, 0, -0.1)),
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "cm"),
              plot.title = element_text(hjust = 0.5, size = 5),
              plot.margin = margin(0,0,0,-10)
              #plot.title = element_blank()
              ) +
        labs(y = "Mean expression", title = paste0("Projection on ", header_name))


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
                legend.margin = margin(0, 0, 0, -8),
                legend.text = element_text(size = 4, margin = margin(0, 0, 0, -0.1)),
                legend.key.size = unit(0.2, "cm"),
                legend.position = "right",
                legend.title = element_blank(),
                plot.title = element_blank()) + 
            labs(y = paste0("Mean % change\nin expr\n(", names(sea_means)[1], " minus ", names(sea_means)[2], ")"))

    # plots[[6]] <- mapply(\(x, y){
    #   highlight_cts <- subdiffallabs |>
    #     dplyr::filter(highlight,
    #                   modno == inds[i]) |>
    #     dplyr::pull(subclass)
    #   df <- x[rownames(x) %in% mod_bc[[inds[i]]], colnames(x) %in% highlight_cts, drop = F] 
    #   out <- data.frame("type" = y, df) |>
    #     pivot_longer(!type, names_to = "ct", values_to = "vals")
    #   return(out)
    # }, sea_means[1:2], names(sea_means)[1:2], SIMPLIFY = F) |>
    #   do.call(what = "rbind") |>
    #   ggpaired(x = "type", 
    #            y = "vals", 
    #            color = "type", 
    #            line.color = "gray", 
    #            line.size = 0.4,
    #            point.size = 0.4,
    #            notch = T) +
    #     facet_wrap(~ct, nrow = 1) + 
    #     labs(x = "", y = "Mean expression") +
    #     theme(legend.position = "none", 
    #         axis.text.x = element_text(size = 4),
    #         axis.title.x = element_blank(),
    #         axis.title.y = element_text(size = 4),
    #         axis.text.y = element_text(size = 4),
    #         strip.text = element_text(size = 4, margin = margin(-5, 0, -5, 0)))  
    return(plots)
  })

  # Plot individual snapshots
  #save_indiv <- file.path(save_path, "indiv")
  #if(!dir.exists(save_indiv)){dir.create(save_indiv, recursive = T)}
  outalllist <- mapply(\(x, y){
    # p <- cowplot::plot_grid(plotlist = x, ncol = 6, align = "h", axis = "b", rel_widths = c(0.6, 0.7, 0.8, 0.8, 0.8, 0.4)) +
    #   coord_cartesian(clip = "off")
    p1 <- cowplot::plot_grid(plotlist = x[1:2], ncol = 3, align = "h", axis = "b", rel_widths = c(0.6, 0.7)) +
      coord_cartesian(clip = "off")
    p2 <- cowplot::plot_grid(plotlist = x[3:5], ncol = 3, align = "h", axis = "b", rel_widths = c(0.6, 0.7, 0.8)) +
      coord_cartesian(clip = "off")
    return(list(p1, p2))
    #p <- cowplot::plot_grid(p1, p2, nrow = 2) 
    #p <- cowplot::plot_grid(plotlist = x, ncol = 3, nrow = 2, align = "h", axis = "b")
    #ggsave(p, file = file.path(save_indiv, paste0(sprintf("%04d", as.numeric(y)), ".pdf")),  height = 2, width = 7, bg = "white", limitsize = F)
    #ggsave(p, file = file.path(save_indiv, paste0(sprintf("%04d", as.numeric(y)), ".svg")),  height = 2, width = 7, bg = "white", limitsize = F)
   # ggsave(p, file = file.path(save_indiv, paste0(y, ".png")),  height = 2, width = 9, bg = "white", limitsize = F)
  }, out_plot, inds, SIMPLIFY = F)
  return(outalllist)

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
                           comp_names = NULL,
                           header_names
                           ){
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = T)
  }

  allpos <- list()
  allneg <- list()
  for(copa_dir in seq_along(copa_dir1)){
    # Load objects
    euc_dist <- list("all"=qread(paste0(copa_dir1[[copa_dir]], "/euclidean_distances/euclidean_sigmods_all.qs")),
                    "pos"=qread(paste0(copa_dir1[[copa_dir]], "/euclidean_distances/euclidean_sigmods_positive.qs")),
                    "neg"=qread(paste0(copa_dir1[[copa_dir]], "/euclidean_distances/euclidean_sigmods_negative.qs")))
    if(is.null(sn_summary_path)){
      sea_means <- qread(paste0(copa_dir1[[copa_dir]], "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
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
    seadist <- fread(data.table = F, file = paste0(copa_dir1[[copa_dir]], "/euclidean_distances/euclidean_distances_all_Subclass.csv"))
    rownames(seadist) <- these_mods
    all_plots <- list(qread(file.path(copa_dir1[[copa_dir]], "sn_proj_objects", "expr_object_bc.qs")), 
                      qread(file.path(copa_dir1[[copa_dir]], "sn_proj_objects", "gsea_object.qs")))

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
        dplyr::filter(mod %in% subset$pos)
      copacohneg_complete <- copacohneg_complete |>
        dplyr::filter(mod %in% subset$neg)
    }

    # Gather percentage differences (pos)
    if(length(copacohpos_complete$mod) > 0){
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
    } else {
      sub_diff_mod <- data.frame()
    }

    # Gather percentage differences (neg)
    if(length(copacohneg_complete$mod) > 0){
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
    } else {
      sub_diff_modneg <- data.frame()
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
        allpos[[copa_dir]] <- plotmods(subdiffallabs = sub_diff_mod1,
                subdiffallper = sub_diff_mod |> mutate(meanpd = meanpd * 100), 
                adconmeans = adconmeans,
                indsmat = copacohpos_complete,
                save_path = file.path(save_path, "pos"),
                file_name = paste0("modList_higherIn", names(sea_means)[1], "_filtered"),
                all_plots = all_plots,
                sea_means = sea_means,
                header_name = header_names[[copa_dir]]
                )
      }
        
      if(length(copacohneg_complete$mod) > 0){
        allneg[[copa_dir]] <- plotmods(subdiffallabs = sub_diff_modneg1,
                  subdiffallper = sub_diff_modneg |> mutate(meanpd = meanpd * 100), 
                  adconmeans = adconmeans,
                  indsmat = copacohneg_complete,
                  save_path = file.path(save_path, "neg"),
                  file_name = paste0("modList_higherIn", names(sea_means)[2], "_filtered"),
                  all_plots = all_plots,
                  sea_means = sea_means,
                  header_name = header_names[[copa_dir]]
                  )
      }
    }
  }
  
  save_indiv <- file.path(save_path, "indiv")
  # Save individual snapshots for pos
  if(length(allpos) > 0){
    for(i in seq_along(allpos[[1]])){
      if(length(allpos) > 1){
        psubl <- lapply(allpos, \(x) x[[i]][[2]])
        psub <- plot_grid(plotlist = psubl, nrow = length(allpos), labels = 'AUTO', label_size = 10)
      } else {
        psub <- allpos[[1]][[i]][[2]]
      }
      rh <- length(allpos) + 1
      p <- cowplot::plot_grid(allpos[[1]][[i]][[1]], psub, nrow = 2, rel_heights = c(1/rh, (rh-1)/rh))
      #ggsave(p, file = file.path(save_indiv, paste0(sprintf("%04d", as.numeric(y)), ".pdf")),  height = 2, width = 7, bg = "white", limitsize = F)
      #ggsave(p, file = "~/test/test.pdf", height = 3, width = 8)
      if(!dir.exists(file.path(save_indiv, "pos")))
        dir.create(file.path(save_indiv, "pos"), recursive = T)
      ggsave(p, file = file.path(save_indiv, "pos", paste0(sprintf("%04d", as.numeric(i)), ".pdf")),  height = length(allpos) + 1, width = 8, bg = "white", limitsize = F)
      ggsave(p, file = file.path(save_indiv, "pos", paste0(sprintf("%04d", as.numeric(i)), ".png")),  height = length(allpos) + 1, width = 8, bg = "white", limitsize = F)
    }
  }

  # Save individual snapshots for neg
  if(length(allneg) > 0){
    for(i in seq_along(allneg[[1]])){
      if(length(allneg) > 1){
        psubl <- lapply(allneg, \(x) x[[i]][[2]])
        psub <- plot_grid(plotlist = psubl, nrow = length(allneg), labels = 'AUTO', label_size = 10)
      } else {
        psub <- allneg[[1]][[i]][[2]]
      }
      rh <- length(allneg) + 1
      p <- cowplot::plot_grid(allneg[[1]][[i]][[1]], psub, nrow = 2, rel_heights = c(1/rh, (rh-1)/rh))
      #ggsave(p, file = file.path(save_indiv, paste0(sprintf("%04d", as.numeric(y)), ".pdf")),  height = 2, width = 7, bg = "white", limitsize = F)
      #ggsave(p, file = "~/test/test.pdf", height = 3, width = 8)
      if(!dir.exists(file.path(save_indiv, "neg")))
        dir.create(file.path(save_indiv, "neg"), recursive = T)
      ggsave(p, file = file.path(save_indiv, "neg", paste0(sprintf("%04d", as.numeric(i)), ".pdf")),  height = length(allneg) + 1, width = 8, bg = "white", limitsize = F)
      ggsave(p, file = file.path(save_indiv, "neg", paste0(sprintf("%04d", as.numeric(i)), ".png")),  height = length(allneg) + 1, width = 8, bg = "white", limitsize = F)
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

# Find overlapping dCoPA mods between two datasets
find_overlapping_mods <- function(paths,
                                  pathnames,
                                  common_pool = NULL,
                                  save_path){
  # paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized"),
  #           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_PFC"))
  # pathnames <- c("bulk_megaset", "bulk_MIT")
  
  # Load lists of significant modules by subclass
  countdf_unfilt <- lapply(paths, \(d){
    compname <- list.files(d)
    modlist <- lapply(compname, \(x){
      path <- file.path(d, x, "modlist_filtered.qs")
      if(file.exists(path)){
        df <- qread(file.path(d, x, "modlist_filtered.qs"))
        outlist <- mapply(\(a, b){
          if(nrow(a) > 0){
            out <- a[, 1:2]
            out$type = b
            return(out)
          } else {
            return(data.frame())
          }
        }, df, c("pos", "neg"), SIMPLIFY = F)
        outlist <- do.call(rbind, outlist)
      } else {
        outlist <- data.frame()
      }
      return(outlist)
    })
    names(modlist) <- compname
    return(modlist)
  })
  names(countdf_unfilt) <- pathnames

  modpersub <- lapply(countdf_unfilt, \(x){
    out <- lapply(x, \(y){
      unique(y$mod)
    })
  })

  # Find intersecting modules (regardless of CT or direction)
  modpersub <- mapply(\(x, y){
    intersect(x, y)
  }, modpersub[[1]], modpersub[[2]], SIMPLIFY = F)

  countdf <- lapply(countdf_unfilt, \(x){
    out <- mapply(\(x, y, z){
      out1 <- x |> 
        mutate("comp" = y) |>
        dplyr::filter(mod %in% z)
      return(out1)
    }, x, 
      names(x), modpersub, SIMPLIFY = F)
    return(out)
  }) 

  # Filter for matching direction
  countdfjoin <- mapply(\(x, y){
    inner_join(x, y, by = join_by("mod" == "mod"), relationship = 'many-to-many') |>
      dplyr::filter(type.x == type.y)
  }, countdf[[1]], countdf[[2]], SIMPLIFY = F)

  ## total # of mods that are:
  # - shared between SEA and MIT
  # - match direction (pos/neg)
  # unlist(lapply(countdfjoin, \(x) length(unique(x$mod))))
    #  conVAll  conVEarly earlyVLate 
    #      116         44         17 

  # How many celltypes on average are significant among the above mods?
  # cts_per_mod <- lapply(countdfjoin, \(z){
  #   temp <- tapply(z$subclass.x, z$mod, list)
  #   return(lapply(temp, \(x) unique(x)))
  # })
  # lapply(cts_per_mod, \(x){
  #   lapply(x, length) |> unlist() |> summary()
  # })
  # $conVAll
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   1.000   2.000   5.000   5.019   7.000  15.000 

  # $conVEarly
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   1.000   1.000   2.000   3.121   4.000  15.000 

  # $earlyVLate
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   1.000   1.000   1.000   1.429   1.500   3.000 

  # Which celltypes are implicated?
  # mods_per_ct <- lapply(countdfjoin, \(z){
  #   zpos <- z |> dplyr::filter(type.x=="pos")
  #   zneg <- z |> dplyr::filter(type.x=="neg")

  #   temppos <- tapply(zpos$mod, zpos$subclass.x,  list)
  #   tempneg <- tapply(zneg$mod, zneg$subclass.x,  list)
  #   return(list("pos" = lapply(temppos, unique),
  #               "neg" = lapply(tempneg, unique)))
  # })
  # lapply(mods_per_ct, \(x){
  #   lapply(x, \(y) lapply(y, length) |> unlist() |> sort(decreasing = T))
  # })

  # $conVAll        
  # $conVAll$pos    
  # Endothelial 
  #           2                                                                                                                                                            
  # $conVAll$neg                                                                                                                                               
  #           L4 IT         L5/6 NP           Pvalb             Sst           L5 IT 
  #             108              96              93              75              72 
  #      L6 IT Car3           L5 ET           Lamp5            Pax6      Lamp5 Lhx6 
  #              65              61              60              34              25 
  #             Vip           L6 IT            Sncg           L6 CT         L2/3 IT 
  #              22              17              14              13              12 
  #             L6b             OPC       Astrocyte       Sst Chodl     Endothelial 
  #               8               6               4               4               1 
  # Oligodendrocyte 
  #               1   
  # $conVEarly               
  # $conVEarly$pos
  #           L4 IT           Pvalb             Sst         L5/6 NP           L5 ET 
  #              75              41              41              39              30 
  #           L5 IT      L6 IT Car3           Lamp5            Pax6           L6 CT 
  #              25              24              18              10               9 
  #             Vip           L6 IT             OPC            Sncg       Sst Chodl 
  #               8               5               4               4               4 
  #       Astrocyte         L2/3 IT            VLMC             L6b      Lamp5 Lhx6 
  #               3               3               3               2               1 
  # Oligodendrocyte 
  #               1 
  # $conVEarly$neg
  #      Chandelier           L6 CT           L6 IT     Endothelial         L2/3 IT 
  #               2               2               2               1               1 
  #             L6b      Lamp5 Lhx6 Oligodendrocyte            Sncg 
  #               1               1               1               1 
  # $earlyVLate
  # $earlyVLate$pos
  #       L6 CT        VLMC   Astrocyte Endothelial         OPC       Pvalb 
  #           2           2           1           1           1           1 

  # $earlyVLate$neg
  #  L6b Pax6 
  #    1    1 


  # #sea celltypes:
  # lapply(countdfjoin, \(x) x$subclass.x) |> unlist() |> unique()
  #  [1] "Endothelial" "L4 IT"       "L5 ET"       "L5 IT"       "L5/6 NP"    
  #  [6] "Lamp5"       "Pvalb"       "Sst"         "L6 IT"       "L6 IT Car3" 
  # [11] "L6 CT"       "Pax6"        "Astrocyte"   "OPC"         "Lamp5 Lhx6" 
  # [16] "Sncg"        "Vip"         "L6b"         "Sst Chodl"   "L2/3 IT"    
  # [21] "VLMC"        "Chandelier" 
  # # mit celltypes:
  # lapply(countdfjoin, \(x) x$subclass.y) |> unlist() |> unique()
  #  [1] "SMC"              "End"              "Inh VIP"          "Exc L4-5 IT-2"   
  #  [5] "Exc L6 IT"        "Exc L3-4 IT"      "Exc L4-5 IT-1"    "Exc L5/6 IT Car3"
  #  [9] "Exc L5 ET"        "Inh PVALB"        "Exc L3-5 IT"      "Inh LAMP5"       
  # [13] "Exc L6b"          "Exc L2-3 IT"      "OPC"              "Inh PAX6"        
  # [17] "Exc L5/6 NP"      "Ast"              "Exc EC"           "Oli"             
  # [21] "Exc L5-6 IT"      "Inh SST"          "Exc L6 CT"        "Per" 
  # common_pool <- list(
  #   c("Endothelial", "SMC", "VLMC", "End", "Per"),
  #   c("L4 IT", "Exc L4-5 IT-2", "Exc L3-4 IT","Exc L4-5 IT-1"),
  #   c("L5 ET",  "Exc L5 ET"),
  #   c("L5 IT", "Exc L4-5 IT-2", "Exc L4-5 IT-1","Exc L3-5 IT", "Exc L5-6 IT"),
  #   c("L5/6 NP", "Exc L5/6 NP"),
  #   c("Lamp5", "Inh LAMP5"),
  #   c("Pvalb", "Inh PVALB"),
  #   c("Sst", "Inh SST"),
  #   c("L6 IT", "Exc L5-6 IT"),
  #   c("L6 IT Car3", "Exc L5/6 IT Car3"),
  #   c("L6 CT", "Exc L6 CT"),
  #   c("Pax6","Inh PAX6"),
  #   c("Astrocyte", "Ast"),
  #   c("OPC", "OPC"),
  #   c("Vip",  "Inh VIP"),
  #   c("L6b", "Exc L6b"),
  #   c("L2/3 IT", "Exc L2-3 IT")
  # )
  
  # Filter to mods where subclasses match
  if(is.null(common_pool)){
    common_pool <- lapply(countdfjoin, \(x){
      c(x$subclass.x, x$subclass.y) |> unique()
    }) |> unlist()
  }

  countdfjoinfilt <- lapply(countdfjoin, \(x){
    these <- lapply(common_pool, \(y){
      which(x$subclass.x %in% y & x$subclass.y %in% y)
    }) |> unlist()
    return(x[these, ])
  })

#  lapply(countdfjoinfilt, \(x) length(unique(x$mod)))
  # $conVAll
  # [1] 89

  # $conVEarly
  # [1] 17

  # $earlyVLate
  # [1] 8

  qsave(countdfjoinfilt, file = save_path)


  # Calculate overlap p-val using Fisher's exact test
  fisherTest_modoverlap <- function(dat1, # SEAAD2024
                                    dat2, # MIT
                                    shared.in.mod, # intersection between SEAAD2024 and MIT
                                    all = 1:1023){ # Fisher's test function
    #total.shared = length(intersect(all,dat1))
    #shared.in.mod = length(intersect(dat2,dat1))
    shared.out.mod = length(dat1) - shared.in.mod
    in.mod.not.shared = length(dat2) - shared.in.mod
    out.mod.not.shared = length(all) - length(dat2) - shared.out.mod
    fisher.test(matrix(c(shared.in.mod,
                        in.mod.not.shared,
                        shared.out.mod,
                        out.mod.not.shared), ncol=2),alternative="greater")$p.val
  }

  fisher_pvals <- lapply(seq_along(countdf_unfilt[[1]]), \(i){
    fisherTest_modoverlap(unique(countdf_unfilt[[1]][[i]]$mod), # SEAAD2024
                          unique(countdf_unfilt[[2]][[i]]$mod), # MIT
                          length(unique(countdfjoinfilt[[i]]$mod)), # intersection between SEAAD2024 and MIT
                          all = 1:1023)
  })
  names(fisher_pvals) <- names(countdf_unfilt[[1]])
  return(fisher_pvals)

  # # conVAll
  
  # # [1] 1.843936e-19

  # # conVEarly
  # fisherTest_modoverlap(unique(countdf_unfilt$bulk_megaset$conVEarly$mod), # SEAAD2024
  #                       unique(countdf_unfilt$bulk_MIT$conVEarly$mod), # MIT
  #                       length(unique(countdfjoinfilt$conVEarly$mod)), # intersection between SEAAD2024 and MIT
  #                       all = 1:1023)
  # # [1] 0.9897561

  # # earlyVLate
  # fisherTest_modoverlap(unique(countdf_unfilt$bulk_megaset$earlyVLate$mod), # SEAAD2024
  #                       unique(countdf_unfilt$bulk_MIT$earlyVLate$mod), # MIT
  #                       length(unique(countdfjoinfilt$earlyVLate$mod)), # intersection between SEAAD2024 and MIT
  #                       all = 1:1023)
  # # [1] 0.5459825


}

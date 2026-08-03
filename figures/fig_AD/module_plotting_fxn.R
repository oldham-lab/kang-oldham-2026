# From ad_associated_gene_module_plotting.R

library(data.table)
library(tidyverse)
library(qs)
library(cowplot)
library(future.apply)
library(ggpubr)
plan(multisession, workers=10)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "generic_enrichment_function.r"))
options(bitmapType = 'cairo')

plot_mods <- function(ind,
                      save_path,
                      file_name){

    save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
    module_output_dir <- save_dir1

  all_plots = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_bc.qs")), 
                qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs")))

  # for indexing and file names
  datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
  if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
  mods <- tapply(datkme[,2], datkme[,3], list)
  modulelengths <- unlist(lapply(mods,length))
  filter_under <- 3
  these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])

  cols <- RColorBrewer::brewer.pal(6, "Paired")
  cols2 <- RColorBrewer::brewer.pal(10, "Spectral")

  out_plot <- lapply(seq_along(ind), \(i){
    j <- which(these_mods == ind[i])
    newtheme <- theme_light() + theme(axis.text.x = element_text(size = 6), 
                                   axis.text.y = element_text(size = 6), 
                                   axis.title.x = element_text(size = 6), 
                                   axis.title.y = element_text(size = 6), 
                                   legend.text = element_text(size = 6),
                                   plot.title = element_blank(),
                                   plot.subtitle = element_blank(),
                                   legend.key.size = unit(0.2, "cm"), 
                                   legend.title = element_blank())
    plots <- lapply(all_plots, \(x) x[[j]] + newtheme)

    plots[[1]]$layers[[1]] <- NULL
    plots[[1]] <- plots[[1]] + 
        geom_line(linewidth = 0.2) +
        labs(x = "Sample") + 
        theme(axis.title.x = element_text(size = 6), 
              axis.title.y = element_text(size = 6),
              legend.text = element_text(size = 7), 
              legend.box.margin = margin(0, 0, 0, -10)) +
        scale_color_manual(values = cols2)
        

    plots[[2]] <- all_plots[[2]][[j]]$data |> 
        dplyr::filter(pval>0) |>
        dplyr::mutate(SetName=factor(SetName, levels=rev(unique(SetName)))) |>
        ggplot(aes(x = SetName, y = pval)) +
            theme_light() +
            geom_bar(stat="identity") +
            theme(axis.text.x = element_text(size=6,hjust=1,vjust = 0.5),
                axis.title.x = element_text(size=6),
                axis.title.y = element_text(size=6),
                axis.text.y = element_text(size=6)) +
            labs(x="", y=bquote(-log[10](p-val))) #+
            #scale_x_discrete(limits = rev(levels(SetName))) +
            #geom_hline(yintercept = gsea_cutFDR, color = "red")
    plots[[2]]$layers[[2]] <- all_plots[[2]][[j]]$layers[[2]]
    plots[[2]] <- plots[[2]] + coord_flip()
    #p <- plot_grid(plotlist = plots, ncol = 4, nrow = 1, align = "h", axis = "bt", rel_widths=c(0.5, 0.85, 0.75, 1, 0.8))
    #return(p)
    return(plots)
  })

    out_plot2 <- lapply(c(1,2), \(x){
        plist <- lapply(out_plot, \(y) y[[x]])
        return(plot_grid(plotlist = plist, nrow = length(plist), align = "v", axis = "rl"))
    })

    outall2 <- plot_grid(plotlist = out_plot2, ncol = 3, align = "h", axis = "bt", rel_widths = c(1,1,0.8))
    ggsave(outall2, file = file.path(save_path, file_name),  height = length(ind) * 1.6, width = 10, bg = "white")
}
# Gather lists of ad-associated genes and examine enrichment in modules
# /home/gugene/data_other/ad_genesets/

library(data.table)
library(tidyverse)
library(qs)
library(cowplot)
library(future.apply)
library(ggpubr)
plan(multisession, workers=10)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "generic_enrichment_function.r"))
options(bitmapType = 'cairo')

# Load modules (bulk)
bulk_seed <- qs::qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
bulk_bc <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
    (\(x) split(x$Gene, x$topmodposbc))()
sn_genelist <- rownames(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
all_genes <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)[,2]

# Load lists of ad risk genes
adlist <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "ad_genesets/five_ad_genesets.qs"))

# Load significant modules (shared between Morabito and SEAAD2024)
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/all_module_index_list.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/neg_module_index_list.qs")))
sigi <- euc_dist$all[[1]] |> unlist() |> unique() # all, subclass, nom - the most enriched mods don't show up in sig mods with fdr cutoff
#sigi <- euc_dist$neg[[1]] |> unlist() |> unique() # all, subclass, nom - but they do show up with nom cutoff
#sigi <- euc_dist$neg[[1]] |> unlist() |> unique() # all, subclass, nom - but they do show up with nom cutoff

# Load significant modules (only SEAAD2024)
euc_dist_sea <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_negative.qs")))
sigi_sea <- euc_dist_sea$all[[1]] |> unlist() |> unique() # all, subclass, nom - the most enriched mods don't show up in sig mods with fdr cutoff
sigi_sea_bc <- euc_dist_sea$all[[2]] |> unlist() |> unique() # all, subclass, fdr

### Selecting modules for plotting
# Run enrichments for all mods 
enrall <- GSHG_custom(bulk_bc, 
                      adlist, 
                      all_genes)
# top 10 enriched mods
top10 <- t(enrall[-1]) |>
  apply(2, \(x) order(x)[1:10]) |> as.data.frame() |> setNames(enrall[,1])
apply(top10, 2, \(x) sum(x %in% sigi_sea))
# [1] 6 9 8 8 8
apply(top10, 2, \(x) sum(x %in% sigi_sea_bc))
# [1] 3 3 1 2 3
# majority of top 10 are significant modules
top10p <- t(enrall[-1]) |>
  apply(2, \(x) x[order(x)[1:10]]) |> as.data.frame() |>
  setNames(enrall[,1])

# Load module snapshots

save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
module_output_dir <- save_dir1

all_plots = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_bc.qs")), 
                qread(file.path(save_dir1,"sn_proj_objects","expr_object_seed.qs")),
                qread(file.path(save_dir1,"sn_proj_objects","bulkcor_object.qs")),
                qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs")),
                qread(file.path(save_dir1,"sn_proj_objects","proj_object.qs")))

# for indexing and file names
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
filter_under <- 3
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])


#apply(top10, 2, \(x) which(x %in% sigi_sea))
# $ADVP
# [1] 1 2 4 5 7 8
# $`Bellenguez et al`
# [1]  1  3  4  5  6  7  8  9 10
# $`Wightman et al`
# [1]  1  2  3  4  6  8  9 10
# $Phenopedia
# [1]  2  3  4  5  6  7  8 10
# $Alzgene
# [1]  1  3  4  5  7  8  9 10
ind <- top10[,1]
#which(ind %in% these_mods)
#[1] 1 2 3 4 5 7 8 9
ind <- ind[ind %in% these_mods]

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
    plots <- lapply(all_plots[c(1,3:4)], \(x) x[[j]] + newtheme)

    plots[[1]]$layers[[1]] <- NULL
    plots[[1]] <- plots[[1]] + 
        geom_line(linewidth = 0.2) +
        labs(x = "Sample") + 
        theme(axis.title.x = element_text(size = 6), 
              axis.title.y = element_text(size = 6),
              legend.text = element_text(size = 7), 
              legend.box.margin = margin(0, 0, 0, -10)) +
        scale_color_manual(values = cols2)
        

    plots[[2]] <- all_plots[[3]][[j]]$data |> 
      ggplot(aes(x = split, y = cors)) + 
        theme_light() +
        geom_violin(aes(fill=type, alpha=type),position=position_dodge(width = 0.8), linewidth=0.1) +
        geom_boxplot(aes(color=type, alpha=type),notch=T,linewidth=0.1,width=0.1,position=position_dodge(width = 0.8), show.legend=F, outlier.shape=NA) +
        geom_hline(color = "black", yintercept=0, alpha=0.6) +
        labs(x="", y="Pearson correlation") +
        scale_x_discrete(labels = scales::label_wrap(10)) +
        theme(legend.title=element_blank(),
            axis.text.x = element_text(size=6),
            axis.title.y = element_text(size=6),
            axis.text.y = element_text(size=6),
            axis.ticks.x=element_blank(),
            legend.position="bottom", 
            legend.direction="vertical", 
            legend.box.margin=margin(-25,0,0,0),
            legend.key.size=unit(0.2, "cm"),
            legend.text=element_text(size=6)) +
        scale_fill_manual(values=c(cols[2], cols[1])) +
        scale_color_manual(values=c("black", "black")) +
        scale_alpha_manual(values=c(1,0.6)) +
        ylim(NA,1.1) +
        stat_compare_means(aes(group = type), label = "p.signif",method.args = list(alternative = "less"), label.y=1, size=2) +
        guides(color = guide_legend(theme(legend.position="none"))) 


    plots[[3]] <- all_plots[[4]][[j]]$data |> 
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
    plots[[3]]$layers[[2]] <- all_plots[[4]][[j]]$layers[[2]]
    plots[[3]] <- plots[[3]] + coord_flip()

    plots[[4]] <- all_plots[[5]][[j]][[1]]$data |> dplyr::filter(indextype=="log_REI") |>
        dplyr::mutate(split=gsub("Subclass", "", split)) |>
        dplyr::mutate(split=factor(split, levels = split |> unique() |> rev())) |>
        ggplot(aes(x = cell_class_match, y = indices, fill = split)) +
          theme_light() +
          geom_col(position=position_dodge(0.5),width=0.5) +
          theme(legend.position="bottom",
                legend.title=element_blank()) + 
          geom_errorbar(aes(ymin = indices-2*mean_se,
                            ymax = indices+2*mean_se),
                        width=0.2,
                        linewidth=0.3,
                        position = position_dodge(0.5)) +
          theme(plot.title = element_text(size=14, hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5,size=12),
                axis.text.x = element_text(size=4, angle=30, vjust=1, hjust=1, margin = margin(t = -2)),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank(), 
                legend.direction = "horizontal",
                legend.text = element_text(size = 6),
                legend.key.size = unit(0.2, "cm"), 
                legend.box.margin = margin(-10, 0, 0, 0)) +
          labs(x = "", y = "") +
          scale_fill_manual(values = cols[c(1,2)])

    #p <- plot_grid(plotlist = plots, ncol = 4, nrow = 1, align = "h", axis = "bt", rel_widths=c(0.5, 0.85, 0.75, 1, 0.8))
    #return(p)
    return(plots)
})

out_plot2 <- lapply(1:4, \(x){
    plist <- lapply(out_plot, \(y) y[[x]])
    return(plot_grid(plotlist = plist, nrow = length(plist), align = "v", axis = "rl"))
})
outall <- plot_grid(plotlist = out_plot2, ncol = 4, align = "h", axis = "bt")
ggsave(outall, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/ad_enrich/advp_modsBC.png"),  height = length(ind) * 1.6, width = 14, bg = "white")
outall2 <- plot_grid(plotlist = out_plot2[c(1,4,3)], ncol = 3, align = "h", axis = "bt", rel_widths = c(1,1,0.8))
ggsave(outall2, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/ad_enrich/advp_modsBC_noBulkCor.png"),  height = length(ind) * 1.6, width = 10, bg = "white")



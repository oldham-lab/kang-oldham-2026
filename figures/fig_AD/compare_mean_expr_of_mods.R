# Compare mean expression of module genes in con vs ad (or early vs late)

library(qs)
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
options(bitmapType = 'cairo')

# Load modules
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(.) split(.$Gene, .$topmodposbc))()

# Load Seattle data
sea_log <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/sn_summary_tables/allmlist_log.qs"))
means_sub <- lapply(sea_log[1:2], rowMeans)

# Load Morabito data
mor_log <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/Morabito_ABIanno/sn_summary_tables/allmlist_log.qs"))
means_sub_mor <- lapply(mor_log[1:2], rowMeans)

# Load significance lists (subclass FDR)
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_negative.qs"))) |>
                 lapply(\(x) x[grep("Subclass",names(x))]) |>
                 lapply(\(x) x[grep("FDR",names(x))] |> unname() |> unlist(recursive=F))
euc_dist_shared <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/all_module_index_list.qs")),
                "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/neg_module_index_list.qs"))) |>
                lapply(\(x) x[grep("Subclass",names(x))]) |>
                lapply(\(x) x[grep("FDR",names(x))] |> unname() |> unlist(recursive=F))

# Plot module means in morabito vs seattle, split by disease state
# seed
plotdfall <- mapply(\(a,b){
    lapply(mod_seed, \(x){
        lapply(a, \(y) mean(y[names(y) %in% x])) |> as.data.frame()
    }) |> do.call(what="rbind") |> setNames(c("AD", "Control")) |>mutate("index"=1:length(mod_seed), "dataset"=b) 
},list(means_sub, means_sub_mor),c("sea", "mor"), SIMPLIFY=F) |> do.call(what="rbind") |>
    mutate(color_all=index %in% unlist(euc_dist_shared$all),
           color_pos=index %in% unlist(euc_dist_shared$pos),
           color_neg=index %in% unlist(euc_dist_shared$neg),
           dataset=factor(dataset),
           AD_minus_con=AD-Control) |>
    mutate(dataset=forcats::fct_recode(dataset,"Morabito"="mor", "SEAAD2024"="sea")) 
    
# Scatter
pall <- plotdfall |>
    ggplot(aes(x=AD, y=Control)) +
        theme_bw() + 
        #geom_point() + 
        geom_abline(linewidth=0.3)+
        geom_point(data=plotdfall[plotdfall$color_pos,],color="red", size=1) +
        geom_point(data=plotdfall[plotdfall$color_neg,],color="green", size=1) +
        theme(text=element_text(size=30),
              plot.margin=margin(1,1,1,1,"cm"),
              legend.position="bottom",
              legend.title=element_blank()) +
        facet_wrap(~dataset, scales="free", ncol=1,nrow=2) +
        labs(x="Log mean expression\n(AD donors)", y="Log mean expression\n(Control donors)") 
ggsave(pall, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_mean_comparison_shared_all.png"),width=6,height=9)


# Split by ct
plotdf_split <- mapply(\(a,b){
    mapply(\(x,w){
        lapply(a, \(y) colMeans(y[rownames(y) %in% x,])) |> do.call(what="rbind") |> as.data.frame() |> 
            tibble::rownames_to_column(var="dx") |> dplyr::mutate("index"=w, "dataset"=b)
    },mod_seed, 1:length(mod_seed), SIMPLIFY=F) |> do.call(what="rbind") 
},list(sea_log[1:2], mor_log[1:2]),c("sea", "mor"), SIMPLIFY=F) 

# scatter plot (facet_wrap)
plotdf_ind <- lapply(2:25, \(x){
    plotdf_split[[1]][,c(1,x,26)] |>
        mutate("subclass"=colnames(plotdf_split[[1]])[x]) |>
        setNames(c("dx","mean_log","index", "subclass")) |>
        mutate(direction = lapply(index, \(a){
                                if(a %in% euc_dist$pos[[x-1]]){
                                    "pos"
                                } else if(a %in% euc_dist$neg[[x-1]]){
                                    "neg"
                                } else {"none"}
                           }) |> unlist()) |>
        tidyr::pivot_wider(names_from=dx, values_from=mean_log)
}) |> do.call(what="rbind")
p <- plotdf_ind |>  dplyr::filter(direction != "none") |> ggplot() +
        theme_test() + 
        #geom_point(data=plotdf[plotdf$color_pos,],aes(x=Alzheimers, y=Control),color="red", size=1) +
        geom_point(aes(x=Alzheimers, y=Control, color=direction), size=0.5) +
        scale_color_manual(values=c("green", "red")) +#, labels=c("Higher in con", "Higher in AD")) +
        geom_abline(linewidth=0.2) +
        theme(text=element_text(size=10),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="none",
            axis.title.x=element_text(size=16, margin=margin(10,0,0,0)),
            axis.title.y=element_text(size=16, margin=margin(0,10,0,0)),
            legend.title=element_blank(),
            plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=10)) +
        labs(x="Log mean expression\n(AD donors)", y="Log mean expression\n(Control donors)") +
        facet_wrap(~subclass, nrow=4, ncol=6)
ggsave(p, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_mean_comparison_seaad2024_by_subclass.png"),width=9,height=6)

# Calculate p-values
# plotdftest <- plotdf |> dplyr::filter(!is.na(direction)) |>
#     dplyr::filter(subclass=="Astrocyte") 
# plotdftest |> group_by(direction) |> 
#     summarise("meanA"=mean(Alzheimers,na.rm=T), "meanB"=mean(Control,na.rm=T),
#               "wilcox_paired_pval"=wilcox.test(Alzheimers,Control, paired=T, alternative="two.sided")$p.val)
# wilcox.test(plotdftest$Alzheimers[1:19],plotdftest$Control[1:19], paired=T, alternative="two.sided")$p.val

plotdfw <- plotdf |> dplyr::group_by(subclass, direction) |>
    dplyr::summarise("wilcox_paired_pval"=wilcox.test(Alzheimers,Control, paired=T, alternative="two.sided")$p.val)

# boxplot (facet_wrap, split by subclass)
plotdfbox <- lapply(2:25, \(x){
    plotdf_split[[1]][,c(1,x,26)] |>
        mutate("subclass"=colnames(plotdf_split[[1]])[x]) |>
        setNames(c("dx","mean_log","index", "subclass")) |>
        mutate(direction = lapply(index, \(a){
                                if(a %in% euc_dist$pos[[x-1]]){
                                    "pos"
                                } else if(a %in% euc_dist$neg[[x-1]]){
                                    "neg"
                                } else {"none"}
                           }) |> unlist())
}) |> do.call(what="rbind")
p <- plotdfbox |>  
    dplyr::filter(!is.na(direction)) |> 
    mutate(dx=factor(dx), direction=factor(direction, levels=c("pos", "neg", "none"))) |>
    mutate(dx=forcats::fct_recode(dx, "AD donors"="Alzheimers", "Con donors"="Control"),
           direction=forcats::fct_recode(direction, "AD-associated"="pos", "Con-associated"="neg", "Not associated"="none")) |>
    ggplot(aes(x=direction, y=mean_log, fill=dx)) +
        theme_test() + 
        geom_boxplot(outlier.size=0.05, size=0.5) +
        scale_fill_manual(values=c("red", "green")) +#, labels=c("Higher in con", "Higher in AD")) +
        theme(text=element_text(size=10),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="bottom",
            axis.title.x=element_blank(),
            axis.title.y=element_text(margin=margin(0,10,0,0)),
            legend.title=element_blank(),
            axis.text.x=element_text(angle=30,hjust=1,vjust=1),
            plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=10)) +
        labs(y="Log mean expression") +
        guides(fill = guide_legend(nrow = 2)) +
        facet_wrap(~subclass, nrow=4, ncol=6)# +
       # ggpubr::stat_compare_means(label="p.signif",label.y=0.9, paired=T)
ggsave(p, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_mean_comparison_seaad2024_by_subclass_boxplot.png"),width=9,height=6)

# Load p-values
filter_under <- 3 
modulelengths <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
      (\(.) tapply(.[,2], .[,3], list))() |>
      lapply(length) |> unlist()
these_mods <- which(modulelengths>filter_under)
seap <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/p_values_all.qs"))

# Plot dist with p-values
plotdf_pval <- lapply(2:25, \(x){
    plotdf_split[[1]][,c(1,x,26)] |>
        filter(index %in% these_mods) |>
        mutate("subclass"=colnames(plotdf_split[[1]])[x]) |>
        setNames(c("dx","mean_log","index", "subclass")) |>
        mutate(direction = lapply(index, \(a){
                                if(a %in% euc_dist$pos[[x-1]]){
                                    "pos"
                                } else if(a %in% euc_dist$neg[[x-1]]){
                                    "neg"
                                } else {"none"}
                           }) |> unlist()) |>
        tidyr::pivot_wider(names_from=dx, values_from=mean_log) %>%
        mutate(diff = Alzheimers-Control,
               pval=seap[[1]][,colnames(seap[[1]])==colnames(plotdf_split[[1]])[x]],
               direction=factor(direction,levels=c("neg", "pos", "none")))
}) |> do.call(what="rbind")
p <- plotdf_pval |> 
    filter(!is.na(Alzheimers)|!is.na(Control)) |>
    ggplot(aes(x=diff, y=pval, color=direction, size=direction, alpha=direction)) +
        theme_test() + 
        geom_point(size=0.5) +
        scale_color_manual(values=c("green", "red", "grey"), labels=c("Higher in con", "Higher in AD", "Not significant")) +
        theme(text=element_text(size=10),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="bottom",
            axis.title.x=element_text(size=16, margin=margin(10,0,0,0)),
            axis.title.y=element_text(size=16, margin=margin(0,10,0,0)),
            legend.title=element_blank(),
            plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=10)) +
        labs(y="p-value", x="Log mean expression difference\n(AD - control)") +
        facet_wrap(~subclass, nrow=4, ncol=6) +
        scale_alpha_manual(values=c(1,1,0.2), labels=c("Higher in con", "Higher in AD", "Not significant")) +
        guides(alpha = guide_legend(legend.position="none")) 

ggsave(p, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_meandiff_vs_pval_seaad2024_by_subclass.png"),width=9,height=6)






######## Split by ct (bc)
plotdf_split <- mapply(\(a,b){
    mapply(\(x,w){
        lapply(a, \(y) if(sum(rownames(y) %in% x)>1){
            colMeans(y[rownames(y) %in% x,])
        } else {
            y[rownames(y) %in% x,]
        } ) |> do.call(what="rbind") |> as.data.frame() |> 
            tibble::rownames_to_column(var="dx") |> dplyr::mutate("index"=w, "dataset"=b)
    },mod_bc, 1:length(mod_bc), SIMPLIFY=F) |> do.call(what="rbind") 
},list(sea_log[1:2], mor_log[1:2]),c("sea", "mor"), SIMPLIFY=F) 

# scatter plot (facet_wrap)
plotdf_ind <- lapply(2:25, \(x){
    plotdf_split[[1]][,c(1,x,26)] |>
        mutate("subclass"=colnames(plotdf_split[[1]])[x]) |>
        setNames(c("dx","mean_log","index", "subclass")) |>
        mutate(direction = lapply(index, \(a){
                                if(a %in% euc_dist$pos[[x-1]]){
                                    "pos"
                                } else if(a %in% euc_dist$neg[[x-1]]){
                                    "neg"
                                } else {"none"}
                           }) |> unlist()) |>
        tidyr::pivot_wider(names_from=dx, values_from=mean_log)
}) |> do.call(what="rbind")
p <- plotdf_ind |>  dplyr::filter(direction != "none") |> ggplot() +
        theme_test() + 
        #geom_point(data=plotdf[plotdf$color_pos,],aes(x=Alzheimers, y=Control),color="red", size=1) +
        geom_point(aes(x=Alzheimers, y=Control, color=direction), size=0.5) +
        scale_color_manual(values=c("green", "red"), labels=c("Higher in con", "Higher in AD")) +
        geom_abline(linewidth=0.2) +
        theme(text=element_text(size=10),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="bottom",
            axis.title.x=element_text(size=16, margin=margin(10,0,0,0)),
            axis.title.y=element_text(size=16, margin=margin(0,10,0,0)),
            legend.title=element_blank(),
            plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=10)) +
        labs(x="Log mean expression\n(AD donors)", y="Log mean expression\n(Control donors)") +
        facet_wrap(~subclass, nrow=4, ncol=6)
ggsave(p, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_mean_comparison_seaad2024_by_subclass_bc.png"),width=9,height=6)

# boxplot (facet_wrap, split by subclass)
plotdfbox <- lapply(2:25, \(x){
    plotdf_split[[1]][,c(1,x,26)] |>
        mutate("subclass"=colnames(plotdf_split[[1]])[x]) |>
        setNames(c("dx","mean_log","index", "subclass")) |>
        mutate(direction = lapply(index, \(a){
                                if(a %in% euc_dist$pos[[x-1]]){
                                    "pos"
                                } else if(a %in% euc_dist$neg[[x-1]]){
                                    "neg"
                                } else {"none"}
                           }) |> unlist())
}) |> do.call(what="rbind")
p <- plotdfbox |>  
    dplyr::filter(!is.na(direction)) |> 
    mutate(dx=factor(dx), direction=factor(direction, levels=c("pos", "neg", "none"))) |>
    mutate(dx=forcats::fct_recode(dx, "AD donors"="Alzheimers", "Con donors"="Control"),
           direction=forcats::fct_recode(direction, "AD-associated"="pos", "Con-associated"="neg", "Not associated"="none")) |>
    ggplot(aes(x=direction, y=mean_log, fill=dx)) +
        theme_test() + 
        geom_boxplot(outlier.size=0.05, size=0.5) +
        scale_fill_manual(values=c("red", "green")) +#, labels=c("Higher in con", "Higher in AD")) +
        theme(text=element_text(size=10),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="bottom",
            axis.title.x=element_blank(),
            axis.title.y=element_text(margin=margin(0,10,0,0)),
            legend.title=element_blank(),
            axis.text.x=element_text(angle=30,hjust=1,vjust=1),
            plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=10)) +
        labs(y="Log mean expression") +
        guides(fill = guide_legend(nrow = 2)) +
        facet_wrap(~subclass, nrow=4, ncol=6)# +
       # ggpubr::stat_compare_means(label="p.signif",label.y=0.9, paired=T)
ggsave(p, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_mean_comparison_seaad2024_by_subclass_boxplot_bc.png"),width=9,height=6)

# Load p-values
filter_under <- 3 
modulelengths <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
      (\(.) tapply(.[,2], .[,3], list))() |>
      lapply(length) |> unlist()
these_mods <- which(modulelengths>filter_under)
seap <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/p_values_all.qs"))

# Plot dist with p-values
plotdf_pval <- lapply(2:25, \(x){
    plotdf_split[[1]][,c(1,x,26)] |>
        filter(index %in% these_mods) |>
        mutate("subclass"=colnames(plotdf_split[[1]])[x]) |>
        setNames(c("dx","mean_log","index", "subclass")) |>
        mutate(direction = lapply(index, \(a){
                                if(a %in% euc_dist$pos[[x-1]]){
                                    "pos"
                                } else if(a %in% euc_dist$neg[[x-1]]){
                                    "neg"
                                } else {"none"}
                           }) |> unlist()) |>
        tidyr::pivot_wider(names_from=dx, values_from=mean_log) %>%
        mutate(diff = Alzheimers-Control,
               pval=seap[[1]][these_mods %in% plotdf_split[[1]]$index,colnames(seap[[1]])==colnames(plotdf_split[[1]])[x]],
               direction=factor(direction,levels=c("neg", "pos", "none")))
}) |> do.call(what="rbind")
p <- plotdf_pval |> 
    filter(!is.na(Alzheimers)|!is.na(Control)) |>
    ggplot(aes(x=diff, y=pval, color=direction, size=direction, alpha=direction)) +
        theme_test() + 
        geom_point(size=0.5) +
        scale_color_manual(values=c("green", "red", "grey"), labels=c("Higher in con", "Higher in AD", "Not significant")) +
        theme(text=element_text(size=10),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="bottom",
            axis.title.x=element_text(size=16, margin=margin(10,0,0,0)),
            axis.title.y=element_text(size=16, margin=margin(0,10,0,0)),
            legend.title=element_blank(),
            plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=10)) +
        labs(y="p-value", x="Log mean expression difference\n(AD - control)") +
        facet_wrap(~subclass, nrow=4, ncol=6) +
        scale_alpha_manual(values=c(1,1,0.2), labels=c("Higher in con", "Higher in AD", "Not significant")) +
        guides(alpha = guide_legend(legend.position="none")) 

ggsave(p, file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_meandiff_vs_pval_seaad2024_by_subclass_bc.png"),width=9,height=6)


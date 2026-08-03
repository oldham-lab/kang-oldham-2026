# Compare number of significant modules before and after filtering SEAAD2024 donors to ones that are similar to Morabito

library(dplyr)
library(ggplot2)
library(qs)
library(data.table)
library(RColorBrewer)

# Load pre-filter mod lists

pre_counts <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/all_module_index_list.qs")),
                   "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                   "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/neg_module_index_list.qs"))) |>
                   lapply(\(x) x[grep("Subclass",names(x))]) |>
                   lapply(\(x) x[grep("FDR",names(x))] |> unname() |> unlist(recursive=F) |> lapply(length) |> unlist())

# Load post-filter mod lists
post_counts <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEA_filtered_to_similar/modules_shared_between_all_datasets/all_module_index_list.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEA_filtered_to_similar/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEA_filtered_to_similar/modules_shared_between_all_datasets/neg_module_index_list.qs"))) |>
                 lapply(\(x) x[grep("Subclass",names(x))]) |>
                 lapply(\(x) x[grep("FDR",names(x))] |> unname() |> unlist(recursive=F) |> lapply(length) |> unlist())

# Create graph objects
p1 <- data.frame("Type"=c("no filtering", "filtered\ndonors"),
                 "Val"=c(pre_counts$all[1], post_counts$all[1])) |>
    mutate(Type=factor(Type, levels=rev(unique(Type)))) |>
    ggplot(aes(y=Type, x=Val,fill=Type)) +
        theme_bw() +
        geom_bar(stat="identity") +
        theme(text=element_text(size=30),
              plot.margin=margin(1,1,1,1,"cm"),
              legend.position="none") +
        labs(y="", x="# of significant mods") +
        scale_fill_manual(values=c("red", "grey"))
ggsave(p1,file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/morabito_vs_sea_before_after_filt_eucDistOverAllCTs.png"), height=3, width=7)


df2 <- list(as.data.frame(pre_counts[2:3]), as.data.frame(post_counts[2:3]))  |> 
            mapply(\(x,y) rownames_to_column(x,"celltype") |> 
                pivot_longer(!celltype,names_to="direction", values_to="count") |>
                mutate(type=y), x=_, c("no filtering", "filtered donors"), SIMPLIFY=F) 

p2 <- do.call(rbind,df2) |>
    mutate(celltype=factor(celltype,levels=rev(unique(celltype))),
           type=factor(type,levels=unique(type)),
           direction=factor(direction, levels=unique(direction))) |>
    mutate(direction=fct_recode(direction,"Higher in AD"="pos", "Higher in con"="neg")) |>
    ggplot(aes(y=celltype, x=count, fill=type)) + 
        theme_bw() +
        geom_bar(stat="identity", position="dodge") +
        theme(text=element_text(size=30),
              plot.margin=margin(1,1,1,1,"cm"),
              legend.position="bottom",
              axis.text.y=element_text(size=15),
              legend.title=element_blank()) +
        facet_wrap(~direction, ncol=2) +
        labs(y="", x="# of significant mods") +
        scale_fill_manual(values=c("grey","red")) 
ggsave(p2,file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/morabito_vs_sea_before_after_filt_by_ct_and_direction.png"), height=8, width=10)

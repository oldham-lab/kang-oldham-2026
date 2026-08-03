# Compare number of significant modules in ROSMAP AD vs CON modules

library(dplyr)
library(ggplot2)
library(qs)
library(data.table)
library(RColorBrewer)

# Load con mod lists
con_counts <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/con_morabito_vs_SEAAD2024/modules_shared_between_all_datasets/all_module_index_list.qs")),
                   "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/con_morabito_vs_SEAAD2024/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                   "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/con_morabito_vs_SEAAD2024/modules_shared_between_all_datasets/neg_module_index_list.qs"))) |>
                   lapply(\(x) x[grep("Subclass",names(x))]) |>
                   lapply(\(x) x[grep("FDR",names(x))] |> unname() |> unlist(recursive=F) |> lapply(length) |> unlist())

# Load AD mod lists
ad_counts <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/all_module_index_list.qs")),
                  "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                  "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/neg_module_index_list.qs"))) |>
                  lapply(\(x) x[grep("Subclass",names(x))]) |>
                  lapply(\(x) x[grep("FDR",names(x))] |> unname() |> unlist(recursive=F) |> lapply(length) |> unlist())


# Create graph objects
p1 <- data.frame("Type"=c("Control", "AD"),
                 "Val"=c(con_counts$all[1], ad_counts$all[1])) |>
    mutate(Type=factor(Type, levels=rev(unique(Type)))) |>
    ggplot(aes(y=Type, x=Val,fill=Type)) +
        theme_bw() +
        geom_bar(stat="identity") +
        theme(text=element_text(size=30),
              plot.margin=margin(1,1,1,1,"cm"),
              legend.position="none") +
        labs(y="", x="# of significant mods") +
        scale_fill_manual(values=c("grey", "red"))
ggsave(p1,file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/rosmap_ad_vs_con_eucDistOverAllCTs.png"), height=3, width=7)


df2 <- list(as.data.frame(con_counts[2:3]), as.data.frame(ad_counts[2:3]))  |> 
            mapply(\(x,y) tibble::rownames_to_column(x,"celltype") |> 
                tidyr::pivot_longer(!celltype,names_to="direction", values_to="count") |>
                mutate(type=y), x=_, c("ROSMAP Con", "ROSMAP AD"), SIMPLIFY=F) 
df2 <- lapply(df2, \(x){
    x$celltype <- factor(x$celltype, levels=x[order(x$count,decreasing=T),]$celltype |> unique() |> rev())
    return(x)
})

p2 <- do.call(rbind,df2) |>
    mutate(type=factor(type,levels=rev(unique(type))),
           direction=factor(direction, levels=unique(direction))) |>
    mutate(direction=forcats::fct_recode(direction,"Higher in AD"="pos", "Higher in con"="neg")) |>
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
        scale_fill_manual(values=c("ROSMAP Con"="grey","ROSMAP AD"="red"), breaks=c("ROSMAP Con", "ROSMAP AD")) +
        guides(fill = guide_legend(nrow = 2)) 

ggsave(p2,file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/rosmap_ad_vs_con_by_ct_and_direction.png"), height=10, width=10)

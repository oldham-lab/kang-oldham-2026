# Function for plotting figure A style figure for copa_compare and copa_conserve outputs
library(data.table)
library(tidyverse)
library(qs)
options(bitmapType = 'cairo')

plot_copa_compare_summary <- function(euc_dist,
                                      plot_title,
                                      pos_neg_vec = c("AD", "con"), # 2 element vector corresponding to "pos" and "neg" (in that order)
                                      save_dir=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/"),
                                      file_suffixes=c("sub_nom", "sub_bc", "sub_fdr", "sup_nom", "sup_bonf", "sup_fdr") # in order of names of euc_dist[[i]]
                                      ){
  if(!dir.exists(save_dir))
    dir.create(save_dir, recursive=T)                                        

  # Collect counts 
  counts_list <- lapply(1:length(file_suffixes), function(t){
    mapply(function(x, y){
      data.frame("ct"=names(x[[t]]), 
                "count"=unlist(lapply(x[[t]],length)),
                "type"=y)
    }, euc_dist, names(euc_dist), SIMPLIFY=F) %>% do.call(rbind,.) %>% 
      dplyr::filter(!(type=="all" & ct!="all")) 
  })

  # Order counts by total significant counts across positive and negative
  counts_order <- lapply(counts_list, function(x){
    x %>% dplyr::filter(type!="all") %>%
      group_by(ct) %>% summarise("sum"=sum(count)) %>% arrange(sum)
  })

  counts_list <- mapply(function(x,y){
    x$ct[x$ct=="all"] <- "Projection over all celltypes"
    x$type[x$type=="all"] <- NA
    x <- x %>%
      mutate(ct=factor(ct, levels=c("Projection over all celltypes",unique(y$ct))))
    return(x)
  }, counts_list, counts_order, SIMPLIFY=F)
 
  # Plot counts
  for(i in seq_along(file_suffixes)){
    p <- ggplot(counts_list[[i]], aes(x=count, y=ct, fill=type))+
      theme_light() +
      labs(x="Number of modules with\nsignificantly different projections",
          y="", title=plot_title) +
      geom_bar(stat="identity") +
      theme(text=element_text(size=30),
            legend.title=element_blank(),
            axis.title.x=element_text(margin=margin(t=20,r=0,b=0,l=0)),
            plot.title=element_text(hjust=0.5)) +
      scale_fill_manual(values = c("#ca5b5b","#50b950"),
                        breaks = c("pos", "neg"),
                        labels = c(paste0("Higher in ", pos_neg_vec[1]),paste0("Higher in ",pos_neg_vec[2])))

    heightvar <- (5/12) * (nrow(counts_list[[i]])-1)/2

    ggsave(p,file=paste0(save_dir,"/panel_A_",file_suffixes[i], ".png"), width=10, height=heightvar, limitsize=F)
    #qsave(p,file=paste0(save_dir,"/panel_A_", file_suffixes[i],".qs"))  
  }

  ## splitbydup
  # Collect counts (subclass, FDR) and distinguish between sig mods that are unique to a celltype and those that are not
  counts_list <- lapply(1:length(file_suffixes), function(t){
    mapply(function(x, y){
      # count the mods that show up in multiple celltypes
      allmods <- unlist(x[[t]])
      dups <- unique(allmods[duplicated(allmods)])

      temp_counts <- lapply(x[[t]], function(x) x[x %in% dups])

      d1 <- data.frame("ct"=names(x[[t]]), 
                "count"=unlist(lapply(x[[t]],function(x) sum(!x %in% dups))),
                "status"="unique",
                "type"=y)
      d2 <- data.frame("ct"=names(x[[t]]), 
                "count"=unlist(lapply(x[[t]],function(x) sum(x %in% dups))),
                "status"="dup",
                "type"=y)
      return(rbind(d1,d2))
    }, euc_dist[2:3], c(paste0("Higher in ", pos_neg_vec[1]), paste0("Higher in ",pos_neg_vec[2])), SIMPLIFY=F) %>% do.call(rbind,.) 
  })

  # Order counts by total significant counts across positive and negative
  order_list <- lapply(counts_list, function(x){
    x %>%
      group_by(ct) %>% summarise("sum"=sum(count)) %>% arrange(sum)
  })

  counts_list <- mapply(function(x,y){
    x %>%
      mutate(ct=factor(ct, levels=c("Projection over all celltypes", unique(y$ct))))
  }, counts_list, order_list, SIMPLIFY=F)

  # Plot counts
  for(i in seq_along(file_suffixes)){
    p <- ggplot(counts_list[[i]], aes(x=count, y=ct, fill=status))+
      theme_light() +
      labs(x="Number of modules with\nsignificantly different projections",
          y="", title=plot_title) +
      geom_bar(stat="identity") +
      theme(text=element_text(size=30),
            #legend.title="Significant in...",
            legend.position="bottom",
            axis.title.x=element_text(margin=margin(t=20,r=0,b=0,l=0)),
            plot.title=element_text(hjust=0.5),
            plot.margin = unit(c(1,1,1,1), "cm")) +
      scale_fill_manual(name="Significant in...",
                        values = c("#ca5b5b","#50b950"),
                        breaks = c("unique", "dup"),
                        labels = c("one celltype", "multiple celltypes")) +
      facet_wrap(~type)

    heightvar <- nrow(counts_list[[i]])/8 + 2

    ggsave(p,file=paste0(save_dir,"/panel_A_splitByDup_", file_suffixes[i], ".png"), width=12, height=heightvar, limitsize=F)  
  }
}



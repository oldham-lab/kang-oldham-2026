library(tidyverse)
library(qs)
library(data.table)
library(scales) ## comma()
library(ggalluvial) ## Ref: https://cheatography.com/seleven/cheat-sheets/ggalluvial/
library(RColorBrewer)
options(bitmapType = 'cairo')

object <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/CHOIR/object_3000varfeatures_prunetree.qs"))
object_meta <- object@meta.data %>% rownames_to_column(var="Cell_ID") %>% dplyr::select(Cell_ID, CHOIR_clusters_0.05)
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE) %>%
  inner_join(object_meta)
inc_mat <- table(cell_annoall$Cluster, cell_annoall$CHOIR_clusters_0.05) %>% as.data.frame


#datinfo <- read.csv("/home/rebecca/SCSN_meta_analysis/datinfo_SCSN_meta_analysis.csv") 
#data_type <- "author_data"
#expr_type <- "counts"
#scaled <- F
#hvgs <- F
#alpha <- .05

#source("/home/rebecca/SCSN_meta_analysis/code/cell_type_stats_fxn.R")

  
#file_path <- paste0("figures/", data_type, "/", expr_type, "/scSHC_test_clusters_riverplot_", alpha, "_FWER_", data_type, "_", expr_type, "_scaled_", scaled, "_HVGs_", hvgs, ".pdf")
file_path <- "~/test/test_river.pdf"  
  
cellinfo <- cell_type_stats(datinfo)
new_clusters <- inc_mat
df <- reshape2::melt(new_clusters, variable.name="New") %>%
    dplyr::filter(value > 0) %>%
    dplyr::select(Var1, Var2, value)
    #dplyr::mutate(New=gsub("new", "", New)) %>%
    #dplyr::arrange(as.numeric(New)) %>%
    #dplyr::mutate(New=paste("New", New)) 

colnames(df) <- c("Original", "New", "No.Nuclei")
  #cellinfo1 <- cellinfo[with(cellinfo, Dataset == names(new_clusters_list)[i]),]
  cellinfo1 <- cell_annoall %>% dplyr::select(Cluster, Cell_Type, Class) %>% 
    dplyr::filter(!duplicated(Cluster))
  df2 <- left_join(df, cellinfo1, by=join_by("Original"=="Cluster"))
  #df <- merge(df, cellinfo1[,c("Cell_Type", "MO_Cell_Class")], 
  #            by.x="Original", by.y="Cell_Type", sort=F)
  #plot_title <- datinfo$Plot_Label[datinfo$Dataset == names(new_clusters_list)[i]]
  #plot_sub <- paste(n_distinct(df$Original), "original clusters,", 
  #                  n_distinct(df$New), "new clusters")
  text_size <- 3
  if(n_distinct(df$Original) > 50){
    text_size <- 2
  }
  # n_row <- 1
  # if(n_distinct(df$MO_Cell_Class) > 6){
  #   n_row <- 2
  # }
for(i in seq_along(unique(df2$Class))){ 
  class <- unique(df2$Class)[i] 
  plotdf <- df2 %>% 
    dplyr::filter(Class==class)
  plot1 <- plotdf %>%
    ggplot(aes(axis1=Original, axis2=New, y=No.Nuclei, fill=Cell_Type)) +
      geom_flow(aes(color=Cell_Type)) + geom_stratum(size=.25) +
      theme_bw() +
      theme(plot.title=element_text(size=16, hjust=0.5),
            plot.subtitle=element_text(size=13, lineheight=1.1, margin=margin(b=7), hjust=0.5),
            legend.title=element_text(size=14),
            legend.text=element_text(size=11),
            legend.position="bottom",
            legend.direction="horizontal",
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.text.x=element_text(size=14, color="black", face="bold", margin=margin(b=5)), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            panel.grid=element_blank(),
            panel.border=element_blank(),
            plot.margin=margin(1, 1, 1, 1, "cm")) +
      geom_text(stat="stratum", aes(label=after_stat(stratum)), size=text_size) +
      labs(title=paste0("Lein DFC ", class),
           subtitle=paste0(length(unique(plotdf$Original))," original clusters, ",length(unique(plotdf$New))," CHOIR clusters")) +
      scale_x_discrete(limits=c("Original", "CHOIR"), expand=c(.1, .1)) +
      scale_y_continuous(trans="sqrt", expand=c(.01, .01)) +
      scale_color_manual(values=brewer.pal(n_distinct(cellinfo1$Cell_Type), "Paired")) +
      scale_fill_manual(values=brewer.pal(n_distinct(cellinfo1$Cell_Type), "Paired")) +
      guides(color="none", fill=guide_legend(title="Subclass", nrow=1, 
                                             override.aes=list(linewidth=.5)))
  ggsave(plot1, file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/CHOIR/"),class,".png"))
}
  

  



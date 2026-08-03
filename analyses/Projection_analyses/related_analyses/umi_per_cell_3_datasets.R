# Graph median UMI per nucleus for Mathys, Morabito, SEAAD2024
# - after subsetting to genes in bulk megaset (~18k genes)

library(ggplot2)
library(RColorBrewer)

brewer_colors <- brewer.pal(6,"Paired")[c(2,4,6)]

umidat <- data.frame("Dataset"=c("Mathys", "Morabito", "SEA-AD 2024"),
                     "median"=c(1474, 6669, 16131))

p <- ggplot(umidat, aes(x=Dataset, y=median)) +
  theme_bw() +
  geom_bar(aes(fill=Dataset),stat="identity") +
  scale_fill_manual(values=brewer_colors) +
  theme(text=element_text(size=30),
        axis.text.x=element_text(angle=30, hjust=1,vjust=1),
        axis.title.y=element_text(size=26),
        plot.margin=margin(1,1,1,1,"cm"),
        legend.position="none") +
  labs(x="", y="Median UMI/nucleus\n(n=18,193 genes)")

ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/median_umi_per_nuc_comparison.png"), width=6)

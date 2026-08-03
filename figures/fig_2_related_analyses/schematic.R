library(ComplexHeatmap)
library(qs)
library(data.table)
library(tidyverse)

# Create randomly generated matrices
rmat <- lapply(1:3, \(x){
  matrix(runif(n = 48, min = 0, max = 10), nrow = 4, ncol = 12)
})

# Create color scales
color_list <- mapply(\(x,y){
  circlize::colorRamp2(c(0, max(x)), c("white", y))
}, rmat, c("red", "green", "blue"), SIMPLIFY = F)

# Draw heatmap borders using global options
ht_opt(
    heatmap_border = TRUE
)
# Create heatmaps
dat_heat <- mapply(\(x, y){
  Heatmap(x,
          col = y,
          cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE
          )
}, rmat, color_list, SIMPLIFY = F)

outdir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/schematic")
if(!dir.exists(outdir))
  dir.create(outdir, recursive = TRUE)

pdf(file.path(outdir, "p1.pdf"))
dat_heat[[1]] %v% dat_heat[[2]] %v% dat_heat[[3]]
dev.off()

# Create pseudobulk mat
pseudo_vec <- lapply(1:5, \(x) sample(1:12, 10))
pseudomat <- lapply(pseudo_vec, \(x){
  df <- do.call(rbind, rmat)
  return(colSums(df[x, ]))
}) 
pseudo_col_ratios <- lapply(pseudo_vec, \(x){
  r <- c(sum(x %in% 1:4)/12, sum(x %in% 5:8)/12, sum(x %in% 9:12)/12)
  return(r/max(r))
})
pseudo_cols <- lapply(pseudo_col_ratios, \(x){
  rgb(x[[1]], x[[2]], x[[3]], max(x))
})
pseudo_scales <- mapply(\(x,y){
  circlize::colorRamp2(c(0, max(x)), c("white", y))
}, pseudomat, pseudo_cols, SIMPLIFY = F)
pseudo_heat <- mapply(\(x, y, z){
  x1 <- matrix(x, nrow = 1)

  Heatmap(x1,
          col = y,
          cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE
          )
}, pseudomat, pseudo_scales, pseudo_col_ratios, SIMPLIFY = F)

pdf(file.path(outdir, "p2.pdf"), width = 10, height = 4)
pseudo_heat[[1]] %v% pseudo_heat[[2]] %v% pseudo_heat[[3]] %v% pseudo_heat[[4]] %v% pseudo_heat[[5]]
dev.off()

# Draw boxplots for color ratios
ratio_df <- mapply(\(x, y){
  data.frame("which" = y, "ratio" = x, "Celltype" = c("A", "B", "C")) |>
    mutate(Celltype = factor(Celltype, levels = unique(Celltype)))
}, pseudo_col_ratios, 1:length(pseudo_col_ratios), SIMPLIFY = F) |>
  do.call("what" = rbind)

p <- ggplot(ratio_df, aes(x = Celltype, y = ratio, fill = Celltype)) +
  theme_classic() + 
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(size = 60),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 24)) +
  scale_fill_manual(values = c("red", "green", "blue")) +
  facet_wrap(~which, ncol = 1)

ggsave(p + theme(legend.position = "none"), 
       file = file.path(outdir, "p3.pdf"), width = 4, height = 7)

leg <- cowplot::get_legend(p)
pdf(file.path(outdir, "legend.pdf"))
grid.newpage()
grid.draw(leg)
dev.off()

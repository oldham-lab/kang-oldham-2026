library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(dendextend)

#####################
# 1st image: representation of bulk megaset using ComplexHeatmap
#####################

# Create ranges for 7 matrices
# range_list <- lapply(1:7, \(x){
#   ind_range <- ((x-1) * 10000):(x * 10000)
#   sample(ind_range, 60000, replace = T) / 10000
# })

# Create dataset size ratios
dat_sizes <- floor(c(255, 182, 274, 145, 190, 63, 597) * 70 / 1706)
dat_names <- c("BrainGVEX", "Brainseq", "CMC", "CMC_HBCC", "GTEx", "NABEC", "ROSMAP")

# Create 7 matrices then concatenate
mat <- lapply(dat_sizes, \(x){
  range <- 1:10000
  input <- sample(range, x * 20)
  matrix(input, nrow = 20, ncol = x)
})
# Create color ranges
range_cols <- RColorBrewer::brewer.pal(7, "Set3")
legend_list <- lapply(range_cols, \(y){
  circlize::colorRamp2(c(1, 10000), c("white", y))
})

# Create heatmaps
h <- mapply(\(x, y){
  Heatmap(x,
          col = y,
          cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE
          )
}, mat, legend_list, SIMPLIFY = F) 
ht_list <- h[[1]] + h[[2]] + h[[3]] + h[[4]] + h[[5]] + h[[6]] + h[[7]]

ht_opt(RESET = TRUE)
ht_opt(
    heatmap_border = FALSE
)
pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_3/schematic/p1.pdf"), width = 10, height = 1.5)
ht_list
#draw(ht_list, ht_gap = unit(0), "mm")
dev.off()

# Create legend
p <- ggplot(data.frame(x = dat_sizes, y = dat_names, fill = dat_names), aes(x = x, y = y, fill = fill)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = range_cols) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(margin = margin(0, 0, 0, -0.05))) 

leg <- cowplot::get_legend(p)
pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_3/schematic/legend.pdf"))
grid.newpage()
grid.draw(leg)
dev.off()

#########
# 2nd image: gene-gene similarity matrix
#########
# Load a simMat
load(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_combined_final.Rdata"))
# Create a simMat
# expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
# expr_t <- t(expr[,-c(1:2)])
# set.seed(20)

# expr_t2 <- cbind(expr_t[,1:10], expr_t[sample(1:nrow(expr_t), nrow(expr_t)), 100:105] - 400)
# input <- cor(expr_t2)

# Plot
# sim_means <- apply(simMat, 2, mean) |> order()
# inds <- c(which(sim_means %in% 1:500),
#           which(sim_means %in% 12000:12500),
#           which(sim_means %in% 9000:9500))

set.seed(26)
inds <- sample(1:ncol(simMat), 20)
input <- simMat[inds, inds]
diag(input) <- 1

pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_3/schematic/p2.pdf"))
Heatmap(input#,
       # show_row_names = FALSE, show_column_names = FALSE
        )
dev.off()
############
# 3rd image: dendrograms
############

# Cluster
rownames(input) <- 1:nrow(input)
colnames(input) <- 1:ncol(input)
cluster1 <- flashClust::flashClust(as.dist(1-input),method="complete")

# Cut
signum <- quantile(input[upper.tri(input)],probs=0.9)
cutree1 <- cutree(cluster1,h = 1-signum)
keptmodsDF <- data.frame(table(cutree1))

# Find initial mods
keptmodscutree <- cutree1[cutree1 %in% keptmodsDF$cutree1[keptmodsDF$Freq >= (3)]]
initialModules <- tapply(as.character(names(keptmodscutree)),factor(keptmodscutree),list)
  
# Plot the cluster
pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_3/schematic/p3.pdf"), height = 3, width = 6)
plot(as.dendrogram(cluster1))
rect.dendrogram(as.dendrogram(cluster1), k=14, lty = 5, lwd = 0, x = c(6,18), col=rgb(0.1, 0.2, 0.4, 0.3), border=0, lower_rect = 0) 
abline(h = 0.31, lty = 2)
dev.off()


## Counts of significant modules in Batiuk 2022

library(data.table)
library(tidyverse)
library(qs)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panel_A_functions.R"))
options(bitmapType = 'cairo')


# Load objects
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/Batiuk_2022/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/Batiuk_2022/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/Batiuk_2022/euclidean_distances/euclidean_sigmods_negative.qs")))

plot_copa_compare_summary(euc_dist,
                          plot_title = "Batiuk 2022",
                          pos_neg_vec = c("SCZ", "con"),
                          save_dir=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/batiuk2022/"))

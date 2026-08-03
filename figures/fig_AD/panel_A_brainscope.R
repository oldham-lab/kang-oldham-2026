## Counts of significant modules in Batiuk 2022

library(data.table)
library(tidyverse)
library(qs)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panel_A_functions.R"))
options(bitmapType = 'cairo')

home <- file.path("/mnt", "bdata", "gugene", "data", "greedy_march_pipeline_output", "finalNonNorm_minsize10_unmerged", "brainSCOPE", "euclidean_distances")

# Load objects
euc_dist <- list("all"=qread(file.path(home,"euclidean_sigmods_all.qs")),
                 "pos"=qread(file.path(home,"euclidean_sigmods_positive.qs")),
                 "neg"=qread(file.path(home,"euclidean_sigmods_negative.qs")))

plot_copa_compare_summary(euc_dist,
                          plot_title = "brainSCOPE",
                          pos_neg_vec = c("SCZ", "con"),
                          save_dir=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/brainSCOPE/"),
                          file_suffixes=c("sub_nom", "sub_bc", "sub_fdr"))

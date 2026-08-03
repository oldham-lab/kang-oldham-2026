# Create a panel A style figure for the copa conserve results 
# (shared modules between mathys, morabito)

library(data.table)
library(tidyverse)
library(qs)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panel_A_functions.R"))
options(bitmapType = 'cairo')

# Load lists of significant modules shared between all three datasets
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/mathys_vs_morabito/modules_shared_between_all_datasets/all_module_index_list.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/mathys_vs_morabito/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/mathys_vs_morabito/modules_shared_between_all_datasets/neg_module_index_list.qs")))

plot_copa_compare_summary(euc_dist,
                          plot_title = "Mathys 2019 vs Morabito",
                          pos_neg_vec = c("AD", "con"),
                          save_dir=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/mathys_vs_morabito/"))


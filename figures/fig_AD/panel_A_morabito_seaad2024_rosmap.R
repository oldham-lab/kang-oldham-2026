# Create a panel A style figure for the copa conserve results 
# (shared modules between the three AD datasets - seaAD, mathys, morabito)

library(data.table)
library(tidyverse)
library(qs)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panel_A_functions.R"))
options(bitmapType = 'cairo')

# Load lists of significant modules shared between all three datasets
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/all_module_index_list.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/neg_module_index_list.qs")))

plot_copa_compare_summary(euc_dist,
                          plot_title = "Morabito vs SEAAD2024 (ROSMAP AD)",
                          pos_neg_vec = c("AD", "con"),
                          save_dir=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/morabito_vs_seaad2024_rosmapAD/"))


# Load lists of significant modules shared between all three datasets
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/con_morabito_vs_SEAAD2024/modules_shared_between_all_datasets/all_module_index_list.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/con_morabito_vs_SEAAD2024/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/COPA_conserve/con_morabito_vs_SEAAD2024/modules_shared_between_all_datasets/neg_module_index_list.qs")))

plot_copa_compare_summary(euc_dist,
                          plot_title = "Morabito vs SEAAD2024 (ROSMAP Con)",
                          pos_neg_vec = c("AD", "con"),
                          save_dir=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/morabito_vs_seaad2024_rosmapCon/"))

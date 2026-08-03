source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/imputation/imputation_fxns.R"))

# scVI
summarise_imputation_results(homedir_impute=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI/"),
                             impute_type="scVI",
                             savedir=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/scVI/")
                             )

# scVI cellbender
summarise_imputation_results(homedir_impute=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI_cellbender/"),
                             impute_type="cb_scVI",
                             savedir=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_scVI/")
                             )

# MAGIC
summarise_imputation_results(homedir_impute=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_MAGIC/"),
                             impute_type="MAGIC",
                             savedir=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/MAGIC/")
                             )

# MAGIC cellbender
summarise_imputation_results(homedir_impute=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_MAGIC_cellbender/"),
                             impute_type="cb_MAGIC",
                             savedir=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_MAGIC/")
                             )

# ALRA
summarise_imputation_results(homedir_impute=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_ALRA/"),
                             impute_type="ALRA",
                             savedir=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/ALRA/")
                             )

# ALRA cellbender
summarise_imputation_results(homedir_impute=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_ALRA_cellbender/"),
                             impute_type="cb_ALRA",
                             savedir=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_ALRA/")
                             )

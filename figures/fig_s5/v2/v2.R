# Previous versions at /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_2_related_analyses/Fig2-lein_pb_indiv_donor_plotting.R

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))

version <- "v2"

########
# Fig. S4: SSv4
########
# Plot DFC
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_SSv4_indiv_donor/")
cell_anno <- fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_ssv4/lein_dfc_metadata.csv"),data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
outdir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s5/"), version, "DFC")
if(!dir.exists(outdir))
  dir.create(outdir, recursive = T)
#donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
donornames <- c("Individual 1", "Individual 2", "Individual 3")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_ssv4/gene_count_means_byDonor.qs"))

plot_fig2_plots_overall(homedir,
                        cell_anno,
                        genemap,
                        outdir,
                        donornames=donornames,
                        san_mean=san_mean,
                        panelc_ymax_r2=1.2,
                        panelc_ymax_rmse=1.2)

# Plot MTG
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_SSv4_indiv_donor_mtg/")
cell_anno <- fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_ssv4/lein_mtg_metadata.csv"), data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
#ctcolvec <- 16 # column of cell_anno that indicates subclass, used for indiv gene graphs for fig1
#sanitygene <- "FBXO43" # gene that corresponds to SANITY threshold in SSv4 data 
outdir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s5/"), version, "MTG")
if(!dir.exists(outdir))
  dir.create(outdir, recursive = T)
#donornames <- c("H200.1023", "H200.1025", "H200.1030") 
donornames <- c("Individual 1", "Individual 2", "Individual 3")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_ssv4/gene_count_means_byDonor.qs"))

plot_fig2_plots_overall(homedir,
                        cell_anno,
                        genemap,
                        outdir,
                        donornames,
                        san_mean=san_mean,
                        panelc_ymax_r2=1.2,
                        panelc_ymax_rmse=1.2)

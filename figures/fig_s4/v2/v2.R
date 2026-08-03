# Previous versions at /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_2_related_analyses/Fig2-lein_pb_indiv_donor_plotting.R

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))

version <- "v2"

###################
# Fig S4
# Run SEAAD2024 DFC con modeling
###################

homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/SEAAD2024con_indiv_donor")
cell_anno <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")) |>
  dplyr::filter(`Overall AD neuropathological Change` == "Not AD")
#sanitygene = "STAMBPL1"
outdir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s4/"), version)
if(!dir.exists(outdir))
  dir.create(outdir, recursive = T)
# donornames <- unique(cell_anno$`Donor ID`)
donornames <- paste0("Individual ", 1:length(unique(cell_anno$`Donor ID`)))
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_con/gene_count_means_byDonor.qs"))

plot_fig2_plots_overall(homedir=homedir,
                        cell_anno=cell_anno,
                        outdir=outdir,
                        donornames=donornames,
                        san_mean=san_mean,
                        panelc_ymax_r2=1.3) # headroom so supertype sanity asterisks aren't clipped

# Only the SVG subfigures feed final figure construction; drop the PNG/PDF
# copies that fxns.R also writes.
file.remove(Sys.glob(file.path(outdir, "part1_*.png")))
file.remove(Sys.glob(file.path(outdir, "part1_*.pdf")))


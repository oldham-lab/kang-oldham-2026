# Previous versions at /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_2_related_analyses/Fig2-lein_pb_indiv_donor_plotting.R

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))

version <- "v3"

########
# Fig. S3: MTG, V1
########
#cell_anno <- fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_10x/lein_mtg_metadata.csv"), data.table=F)
cell_anno <- fread(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_MTG.csv"), data.table = F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/MTG_indiv_donor/")
outdir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s3/"), version, "MTG")
if(!dir.exists(outdir))
  dir.create(outdir, recursive = T)
#donornames <- c("H200.1023", "H200.1025", "H200.1030")
donornames <- c("Individual 1", "Individual 2", "Individual 3")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_10x/gene_count_means_byDonor.qs"))
 
plot_fig2_plots_overall(homedir = homedir,
                        cell_anno = cell_anno,
                        genemap = genemap,
                        outdir = outdir,
                        donornames = donornames,
                        san_mean = san_mean,
                        panelc_ymax_r2 = 1.2)

plot_fig2_plots_indiv(homedir = homedir,
                      cell_anno = cell_anno,
                      genemap=genemap,
                      ctcolvec = 3,
                      outdir=outdir,
                      cellidcol = 2,
                      save_name = "panel_D")

# Run V1
cell_anno <- fread(data.table=F,file=file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_V1.csv"))
cell_anno <- cell_anno[,c(3,1,2,4:8)] # just for consistencygenemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/V1_indiv_donor/")
outdir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s3/"), version, "V1")
if(!dir.exists(outdir))
  dir.create(outdir, recursive = T)
#donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
donornames <- c("Individual 1", "Individual 2", "Individual 3")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_v1_10x/gene_count_means_byDonor.qs"))

plot_fig2_plots_overall(homedir=homedir,
                        cell_anno=cell_anno,
                        genemap=genemap,
                        outdir=outdir,
                        donornames=donornames,
                        san_mean=san_mean)

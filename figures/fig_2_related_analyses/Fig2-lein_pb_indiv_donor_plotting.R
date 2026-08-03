source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/Fig1-functions.R"))

# Run DFC
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/")
cell_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
#sanitygene = "STAMBPL1"
outdir = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/dfc_10x/")
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
san_mean <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data_other/lein_dfc/gene_count_means_byDonor.qs"))
 
plot_fig1_plots_overall(homedir=homedir,
                        cell_anno=cell_anno,
                        genemap=genemap,
                        outdir=outdir,
                        donornames=donornames,
                        san_mean=san_mean,
                        panelc_ymax_r2 = 1.2)
                        
plot_fig1_plots_indiv(homedir=homedir,
                      cell_anno=cell_anno,
                      genemap=genemap,
                      ctcolvec=1,
                      outdir=outdir,
                      cellidcol=3)

# Run MTG
#cell_anno <- fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_10x/lein_mtg_metadata.csv"), data.table=F)
cell_anno <- fread("/mnt/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_MTG.csv", data.table = F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/MTG_indiv_donor/")
outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_2/mtg_10x/")
donornames <- c("H200.1023", "H200.1025", "H200.1030")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_10x/gene_count_means_byDonor.qs"))

plot_fig1_plots_overall(homedir = homedir,
                        cell_anno = cell_anno,
                        genemap = genemap,
                        outdir = outdir,
                        donornames = donornames,
                        san_mean = san_mean,
                        panelc_ymax_r2 = 1.2)

plot_fig1_plots_indiv(homedir = homedir,
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
outdir = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/v1_10x/")
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_v1_10x/gene_count_means_byDonor.qs"))

plot_fig1_plots_overall(homedir=homedir,
                        cell_anno=cell_anno,
                        genemap=genemap,
                        outdir=outdir,
                        donornames=donornames,
                        san_mean=san_mean)

# Plot figure s4
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/")
cell_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x/figs4/")
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means_byDonor.qs"))

inputgenes_unordered <- c("VIP", "SST", "LAMP5", "LHX6", "PAX6", "GAD2", "SLC32A1", "SLC17A6","SLC17A8", "RBFOX3")
inputgenes_ordered <- san_mean[[1]][names(san_mean[[1]]) %in% inputgenes_unordered] %>% sort %>% names

plot_fig1_plots_indiv(homedir=homedir,
                      cell_anno=cell_anno,
                      genemap=genemap,
                      ctcolvec=1,
                      outdir=outdir,
                      cellidcol=3,
                      inputgenevec1 = inputgenes_ordered[1:5],
                      inputgenevec1_names = inputgenes_ordered[1:5],
                      save_name="figs4_part1")

plot_fig1_plots_indiv(homedir=homedir,
                      cell_anno=cell_anno,
                      genemap=genemap,
                      ctcolvec=1,
                      outdir=outdir,
                      cellidcol=3,
                      inputgenevec1 = inputgenes_ordered[6:10],
                      inputgenevec1_names = inputgenes_ordered[6:10],
                      save_name="figs4_part2")

# Run SEAAD2024 con
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/SEAAD2024con_indiv_donor")
cell_anno <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")) |>
  dplyr::filter(`Overall AD neuropathological Change` == "Not AD")
#sanitygene = "STAMBPL1"
outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_2/sea_con/")
donornames <- unique(cell_anno$`Donor ID`)
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_con/gene_count_means_byDonor.qs"))

plot_fig1_plots_overall(homedir=homedir,
                        cell_anno=cell_anno,
                        outdir=outdir,
                        donornames=donornames,
                        san_mean=san_mean)

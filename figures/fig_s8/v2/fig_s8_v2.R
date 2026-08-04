# Figure S8 (v2) panels -- MTG analogue of fig_5/v5 (see assemble_figure_s8_v2.py).
# Identical panel pipeline to fig_5/v5, but on the MTG networks (fig_s8/v1) and MTG REI
# data. Produces:
#   panel_B.svg  - labeled dendrogram + barplots (Fig S8 panel a)
#   panel_C.svg  - Jorstad MTG heatmap, row dendrogram (Fig S8 panel b)
#   panel_D.svg  - Gabitto MTG heatmap, row dendrogram (Fig S8 panel c)
#   + legend, mod_eig.csv, branchpoint tables, panel coords.
# (Panel filenames kept as B/C/D to mirror fig_5; the assembly relabels them a/b/c.)
# The MTG FindModules networks already exist in fig_s8/v1, so nothing is re-run here.

library(ComplexHeatmap); library(tidyverse); library(data.table); library(dendextend); library(qs); library(showtext)
showtext_auto()

save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s8/v2")
if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
# v2 standalone (HGNC re-run): the two MTG FindModules networks are generated in THIS script
# (ported from fig_s8/v1/v1.R) rather than read from a pre-built fig_s8/v1, so v2 no longer
# depends on v1. jor_net / gab_net are resolved by glob after FindModules runs (see below).

# ---- cell-type metadata (MTG) ----
class_info <- fread(data.table=F, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  select(Subclass, Class) |> filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(Class, "Neuronal: GABAergic"~"GABAergic","Neuronal: Glutamatergic"~"Glutamatergic","Non-neuronal and Non-neural"~"Non-neuronal")) |>
  arrange(Subclass) |>
  mutate(Subclass_fixed = factor(c("Astro","Chandelier","Endo","L2/3 IT","L4 IT","L5 ET","L5 IT","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b","Lamp5","Lamp5 Lhx6","Micro/PVM","OPC","Oligo","Pax6","Pvalb","Sncg","Sst","Sst Chodl","VLMC","Vip"),
    levels=c("L2/3 IT","L4 IT","L5 ET","L5 IT","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b","Chandelier","Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Sncg","Sst","Vip","Sst Chodl","Astro","Oligo","OPC","Micro/PVM","Endo","VLMC"))) |>
  arrange(Subclass_fixed)
allcts     <- c("L2/3 IT","L4 IT","L5 ET","L5 IT","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b","Chandelier","Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Astrocyte","Endothelial","Microglia-PVM","Oligodendrocyte","OPC","VLMC")
allcts_cap <- c("L2/3 IT","L4 IT","L5 ET","L5 IT","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b","Chandelier","LAMP5","LAMP5 LHX6","PAX6","PVALB","SNCG","SST","SST CHODL","VIP","Astrocyte","Endothelial","Microglia-PVM","Oligodendrocyte","OPC","VLMC")
class_info[,4] <- allcts_cap[match(class_info[,1], allcts)]

# ---- MTG REI + module filter -> cell-type dendrogram (clustermin) ----
# HGNC fix: Gabitto MTG REI repointed from stale R sn_proj_indices to the HGNC-fixed Python
# SEA-AD output (region MTC = MTG). No `module` column; select(class_info[,1]) keeps the 24 celltypes.
rei1 <- fread(data.table=F, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/mod_means/log_REI/mod_means_Con_bulk_megaset.csv")) |> select(class_info[,1])
rei2 <- fread(data.table=F, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinMTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_con.csv")) |> select(class_info[,3])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table=F, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv")) |> (\(x) tapply(x[,2], x[,3], list))()
these_mods <- which(lapply(mod_bc,length) |> unlist() > 3)
sigcount_bonf <- fread(data.table=F, file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]
mod_seed <- mod_seed[these_mods]
rei1 <- rei1[these_mods, ]; rei2 <- rei2[these_mods, ]; colnames(rei2) <- colnames(rei1)
rownames(rei1) <- 1:nrow(rei1); rownames(rei2) <- 1:nrow(rei2)     # Mod ids 1..1016 (FindModules simMat/expr must match)
# HGNC re-run: recompute zero-variance module positions from the fixed data instead of hardcoding
# c(753,972) — the fix recovered Gabitto's all-histone module 972, so only Jorstad's 753 remains.
zerovar_vec <- sort(union(which(apply(rei1,1,var)==0), which(apply(rei2,1,var)==0)))
rei1 <- rei1[-zerovar_vec, ]; rei2 <- rei2[-zerovar_vec, ]         # zero-variance projections (MTG)
ctm1 <- pmin(cor(rei1), cor(rei2))
clustermin <- hclust(as.dist(1 - ctm1), method="complete") |> as.dendrogram()

# ---- Standalone MTG network generation (ported from fig_s8/v1/v1.R) ----
# mod-mod consensus-min similarity matrix used as the FindModules simMat.
m1 <- cor(t(rei1)); m2 <- cor(t(rei2))
mcm1 <- pmin(m1, m2)
colnames(mcm1) <- paste0("Mod", colnames(mcm1)); rownames(mcm1) <- paste0("Mod", rownames(mcm1))

setwd(file.path(Sys.getenv("FINDMODULES_DIR", "/home/gugene/code/git/FindModules"), "FindModules/R/"))
source("FindModules.R"); source("map_identifiers_function.R"); source("FM_helper_fxns.R")
source("FindModules.R"); source("find_seed_genes_greedy_march_megaset.R"); source("similarityType.R")
source("plotting_functions.R"); source("overlapType.R"); source("networkOutputs.R")
source("module_quant_functions.R"); source("iteration_code.R")
setwd(save_dir)

# Gabitto MTG (HGNC-fixed Python MTC REI)
proj <- fread(data.table=F, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/mod_means/log_REI/mod_means_Con_bulk_megaset.csv")) |> select(class_info[,1])
expr <- data.frame("Gene" = paste0("Mod", 1:length(these_mods)), proj[these_mods, ])
expr <- expr[-zerovar_vec, ]
FindModules(projectname="Gabitto_MTG_consensusMin_noMerge", data_cols=expr, genes=NULL,
  metadata_cols=1, sampleIndex=2:ncol(expr), sampleGroups=NULL, subset=c("None"), simMat=mcm1,
  saveSimMat=FALSE, simType=c("Pearson"), overlapTO=FALSE, TOtype=c("default"), TOdenom=c("default"),
  beta=1, iterate=TRUE, minSizevec=c(3,5,8), greedyMarch=FALSE, signumType=c("rel"),
  signumvec=c(.97,.96,.95), minMEcorvec=c(0.95), merge.by=c("ME"), merge.param=1,
  export.merge.comp=TRUE, loadTree=FALSE, writeKME=TRUE, writeModSnap=TRUE,
  modSnapExprVal=c("meanexpr"), floor="default", prompt_zero_values=F)

# Jorstad MTG (reference; HGNC-invariant)
proj <- fread(data.table=F, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinMTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_con.csv")) |> select(class_info[,3])
expr <- data.frame("Gene" = paste0("Mod", 1:length(these_mods)), proj[these_mods, ])
expr <- expr[-zerovar_vec, ]
FindModules(projectname="Jorstad_MTG_consensusMin_noMerge", data_cols=expr, genes=NULL,
  metadata_cols=1, sampleIndex=2:ncol(expr), sampleGroups=NULL, subset=c("None"), simMat=mcm1,
  saveSimMat=FALSE, simType=c("Pearson"), overlapTO=FALSE, TOtype=c("default"), TOdenom=c("default"),
  beta=1, iterate=TRUE, minSizevec=c(3,5,8), greedyMarch=FALSE, signumType=c("rel"),
  signumvec=c(.97,.96,.95), minMEcorvec=c(0.95), merge.by=c("ME"), merge.param=1,
  export.merge.comp=TRUE, loadTree=FALSE, writeKME=TRUE, writeModSnap=TRUE,
  modSnapExprVal=c("meanexpr"), floor="default", prompt_zero_values=F)

# Resolve the generated networks by glob (lowest signum = largest network, minSize3) — robust to
# the signum shift caused by the HGNC-fixed REI (mirrors fig_5/v5's find_net).
find_net <- function(prefix){
  hits <- Sys.glob(file.path(save_dir, paste0(prefix, "_consensusMin_noMerge_Modules"),
                             "Pearson-no_TO_signum*_minSize3_merge_ME_1_*"))
  if(length(hits) == 0) stop("No minSize3 ", prefix, " network in ", save_dir)
  hits[which.min(as.numeric(sub(".*signum([0-9.]+)_minSize3.*", "\\1", basename(hits))))]
}
jor_net <- find_net("Jorstad_MTG")
gab_net <- find_net("Gabitto_MTG")

# ---- meta-module eigengenes (Jorstad) -> metacluster_min, numbering, barplot counts ----
mod_eig <- fread(data.table=F, file.path(jor_net,"Module_eigengenes.csv")) |> column_to_rownames("Sample")
rownames(mod_eig) <- class_info[,4]; jor_order <- colnames(mod_eig); colnames(mod_eig) <- 1:ncol(mod_eig)
metacluster_min <- hclust(as.dist(1 - cor(mod_eig)), method="complete") |> as.dendrogram()
qsave(metacluster_min, file.path(save_dir,"jorstad_consensusMin_dendro.qs"))
kme <- fread(data.table=F, file.path(jor_net,"kME_table_.csv"))
mod_fdr <- tapply(kme[,1], kme[,6], list) |> lapply(\(x) as.numeric(gsub("Mod","",x)))
mod_fdr_gene <- lapply(mod_fdr, \(x) mod_seed[x] |> unlist() |> unique())
modcountdf  <- data.frame(mod=names(mod_fdr), FDR=lapply(mod_fdr,length)|>unlist()) |>
  pivot_longer(!mod, names_to="sig_cut", values_to="mod_count") |> mutate(mod_per = mod_count/length(these_mods)*100)
genecountdf <- data.frame(mod=names(mod_fdr_gene), FDR=lapply(mod_fdr_gene,length)|>unlist()) |>
  pivot_longer(!mod, names_to="sig_cut", values_to="gene_count") |> mutate(gene_per = gene_count/18913*100)
countdf <- full_join(modcountdf[,-c(2:3)], genecountdf[,-c(2:3)], by=join_by(mod)) |> as.data.frame()
countdf <- countdf[match(jor_order, countdf$mod), ]; countdf$mod <- 1:nrow(countdf)

# ---- Jorstad panel WITHOUT column dendrogram (barplots + heatmap) + geometry for compositing ----
p <- Heatmap(as.matrix(mod_eig), name="Eigenmodule", cluster_rows=clustermin, cluster_columns=metacluster_min,
  show_column_dend=FALSE,
  top_annotation=HeatmapAnnotation("% of all mods"=anno_barplot(countdf[,2], axis_param=list(at=c(0,2))),
                                   "% of all genes"=anno_barplot(countdf[,3], axis_param=list(at=c(0,1)))),
  heatmap_legend_param=list(title_gp=gpar(fontface="plain"), title_position="topcenter"),
  column_title_side="bottom", column_title="Meta-modules", row_title="1 - cor", row_names_side="left",
  column_names_gp=gpar(fontsize=9), row_names_gp=gpar(fontsize=12), show_heatmap_legend=FALSE)
svg(file.path(save_dir,"/panel_C_consensusMin_Jorstad_noDend.svg"), width=15, height=7)
draw(p, padding=unit(c(6,6,6,24),"mm"), heatmap_legend_side="bottom")
decorate_heatmap_body("Eigenmodule", {
  tl <- deviceLoc(unit(0,"npc"),unit(1,"npc")); tr <- deviceLoc(unit(1,"npc"),unit(1,"npc"))
  x0 <<- convertX(tl$x,"inches",valueOnly=TRUE); x1 <<- convertX(tr$x,"inches",valueOnly=TRUE)
  ybody <<- convertY(tl$y,"inches",valueOnly=TRUE) })
decorate_annotation("% of all mods", { tl <- deviceLoc(unit(0,"npc"),unit(1,"npc")); ybar <<- convertY(tl$y,"inches",valueOnly=TRUE) })
dev.off()
writeLines(c(paste0("W_px=",15*72), paste0("H_px=",7*72), paste0("x0=",x0*72), paste0("x1=",x1*72),
             paste0("ybar_top=",(7-ybar)*72), paste0("ybody_top=",(7-ybody)*72)),
           file.path(save_dir,"jorstad_panel_coords.txt"))

# colour legend
heatmap_legend <- Legend(col_fun=p@matrix_color_mapping@col_fun,
  title="Eigenmodules\n(PC1 of REIs\nfor merged\nmodules)", title_gp=gpar(fontface="plain"),
  title_position="topcenter", direction="vertical")
svg(file.path(save_dir,"/panel_C_consensusMin_Jorstad_legend.svg"), width=1.5, height=4); grid.newpage(); draw(heatmap_legend); dev.off()

# ---- labeled-dendrogram inputs: mod_eig.csv (colour names, dendrogram-ordered) + branchpoint tables ----
mod_eig_lab <- fread(data.table=F, file.path(jor_net,"Module_eigengenes.csv")) |> column_to_rownames("Sample")
rownames(mod_eig_lab) <- class_info[,1]
mod_eig_lab <- mod_eig_lab[order.dendrogram(clustermin), order.dendrogram(metacluster_min)]
fwrite(mod_eig_lab, file.path(save_dir,"mod_eig.csv"), row.names=TRUE)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s8/build_branchpoint_table.R"))
bp_table <- build_branchpoint_table(mod_eig_file=file.path(save_dir,"mod_eig.csv"), cut_height=0.3,
                                    output_file=file.path(save_dir,"branchpoint_table_modeig.csv"))
bp_table <- bp_table |> mutate(Elements = map(Elements, ~ str_split(.x, ",\\s*")[[1]]))
numfdrgene <- length(unique(unlist(mod_seed))); nummodgene <- length(mod_seed)
gene_per <- sapply(bp_table$Elements, \(x) length(unlist(mod_fdr_gene[names(mod_fdr_gene) %in% x]))/numfdrgene)
bp_table$Pct_of_Total_genes <- signif(gene_per*100, 3)
mod_per <- sapply(bp_table$Elements, \(x) length(unlist(mod_fdr[names(mod_fdr) %in% x]))/nummodgene)
bp_table$Pct_of_Total <- signif(mod_per*100, 3)
bp_table <- bp_table[, c(1:4, 6, 5)]
fwrite(bp_table, file.path(save_dir,"branchpoint_table_modeig_with_genes.csv"))

# ---- render labeled dendrogram + composite onto Jorstad panel -> panel_B.svg + panel_C_consensusMin_Jorstad.svg ----
render_py  <- Sys.getenv("PYTHON_BIN", "/home/gugene/miniconda3/bin/python")
render_scr <- file.path(save_dir,"make_labeled_dendrogram_only.py")
composite  <- file.path(save_dir,"composite_dendrogram_panelC.py")
system2(render_py, c(shQuote(render_scr),
  shQuote(file.path(save_dir,"mod_eig.csv")), shQuote(file.path(save_dir,"branchpoint_table_modeig_with_genes.csv")),
  "-o", shQuote(file.path(save_dir,"jorstad_dendro_labeled.svg")),
  "--cut-height","0.3","--font-size","14","--fig-width","36","--fig-height","8",
  "--coords-out", shQuote(file.path(save_dir,"jorstad_dendro_coords.txt"))))
system2(render_py, c(shQuote(composite), shQuote(save_dir)))
system2("rsvg-convert", c("-f","pdf","-o", shQuote(file.path(save_dir,"panel_C_consensusMin_Jorstad.pdf")),
  shQuote(file.path(save_dir,"panel_C_consensusMin_Jorstad.svg"))))

# ---- heatmap-only panels: panel_C.svg (Jorstad, no x-axis) + panel_D.svg (Gabitto, x-axis) ----
panel_heatmap <- function(me, out, coords_out, xaxis=TRUE){
  hm <- Heatmap(as.matrix(me), name="Eigenmodule", cluster_rows=clustermin, cluster_columns=metacluster_min,
    show_column_dend=FALSE, column_title_side="bottom", column_title=if(xaxis)"Meta-modules" else NULL,
    show_column_names=xaxis, row_title="1 - cor", row_names_side="left",
    column_names_gp=gpar(fontsize=9), row_names_gp=gpar(fontsize=12), show_heatmap_legend=FALSE)
  svg(out, width=15, height=5); draw(hm, padding=unit(c(6,6,1,24),"mm"))
  decorate_heatmap_body("Eigenmodule", {
    tl <- deviceLoc(unit(0,"npc"),unit(1,"npc")); tr <- deviceLoc(unit(1,"npc"),unit(1,"npc"))
    writeLines(c(paste0("W_px=",15*72), paste0("x0=",convertX(tl$x,"inches",valueOnly=TRUE)*72),
                 paste0("x1=",convertX(tr$x,"inches",valueOnly=TRUE)*72)), coords_out) })
  dev.off()
}
jor_me <- fread(data.table=F, file.path(jor_net,"Module_eigengenes.csv")) |> column_to_rownames("Sample")
rownames(jor_me) <- class_info[,4]; colnames(jor_me) <- 1:ncol(jor_me)
panel_heatmap(jor_me, file.path(save_dir,"panel_C.svg"), file.path(save_dir,"panel_C_coords.txt"), xaxis=FALSE)
gab_me <- fread(data.table=F, file.path(gab_net,"Module_eigengenes.csv")) |> column_to_rownames("Sample") |> select(all_of(jor_order))
rownames(gab_me) <- class_info[,4]; colnames(gab_me) <- 1:ncol(gab_me)
panel_heatmap(gab_me, file.path(save_dir,"panel_D.svg"), file.path(save_dir,"panel_D_coords.txt"), xaxis=TRUE)

message("fig_s8 v2 panels written: panel_B.svg (a), panel_C.svg (b), panel_D.svg (c) + legend")

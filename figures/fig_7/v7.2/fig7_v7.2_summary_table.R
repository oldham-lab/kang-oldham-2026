# v7.2
# Finalized summary table
# (Only produces panels e-f)

library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)
library(showtext)
showtext_auto()

version <- "v7.2/"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)

if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.2/ad_db.R"))

####### Start with DFC
g_overlaps <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.2/dfc_overlaps.csv"))

# Run GSEA
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func_optimized.R"))
b_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = T)
g_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = F)
gsea_out <- rbind(b_out, g_out)

# Calculate FDR
gsea_pval_to_fdr <- function(gsea_out) {                              
    meta_cols <- c("SetID", "SetName", "SetSize")                       
    pval_cols <- setdiff(colnames(gsea_out), meta_cols)                 
    fdr_out <- gsea_out                                                 
    fdr_out[, pval_cols] <- lapply(gsea_out[, pval_cols, drop = FALSE], 
                                    p.adjust, method = "BH")            
    fdr_out                                                             
  }                                                                     
gsea_out_fdr <- gsea_pval_to_fdr(gsea_out)


# Significant genesets per celltype (p < GSEA_PTHRESH)
# GSEA_PTHRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   gsea_out$SetName[!is.na(pvals) & pvals < GSEA_PTHRESH]
# }), celltypes)

# ── AD gene database ────────────────────────────────────────────────────────
library(gt)

dfc_hits    <- fread(file.path(save_dir, "ad_db_summary_table_dfc.csv"))
pubmed_hits <- setNames(as.integer(dfc_hits$Pubmed_total_hits), dfc_hits$Gene)

# ── Build table ─────────────────────────────────────────────────
celltypes  <- unique(g_overlaps$Celltype)

# # Significant genesets per celltype (BH-adjusted p < GSEA_FDR_THRESH)
# GSEA_FDR_THRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   fdr   <- p.adjust(pvals, method = "BH")
#   gsea_out$SetName[!is.na(fdr) & fdr < GSEA_FDR_THRESH]
# }), celltypes)
# qsave(gsea_sig, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/gsea_sig.qs"))

CT_RENAME <- c(
  "Lamp5"      = "LAMP5",
  "Lamp5 Lhx6" = "LAMP5 LHX6",
  "Pax6"       = "PAX6",
  "Pvalb"      = "PVALB",
  "Sncg"       = "SNCG",
  "Sst"        = "SST",
  "Vip"        = "VIP",
  "Sst Chodl"  = "SST CHODL"
)

build_row <- function(ct) {
  genes  <- g_overlaps$Gene[g_overlaps$Celltype == ct]
  n_genes <- length(genes)

  # Top 3 GO gene sets for this celltype (columns assumed in same order as celltypes)
  ct_idx  <- which(celltypes == ct)
  pvals   <- gsea_out[[3 + ct_idx]]
  go_mask <- grepl("^GO", gsea_out$SetName)
  if (!is.null(pvals)) {
    go_pvals <- ifelse(go_mask, pvals, NA)
    sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
    n_sig    <- sum(sig_mask)
    if (n_sig == 0) {
      top_str <- "—"
    } else {
      top_idx  <- order(go_pvals)[1:min(3, n_sig)]
      top_sets <- gsea_out$SetName[top_idx]
      top_pval <- pvals[top_idx]
      top_str  <- paste0(
        seq_along(top_sets), ". ", top_sets,
        " (p=", formatC(top_pval, format = "e", digits = 2), ")",
        collapse = "<br>"
      )
    }
  }

  # AD genes: top 10 by PubMed hits
  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(ad_db_dfc)]
  hits_up     <- genes_upper[genes_upper %in% names(ad_db_dfc)]

  if (length(hits) > 0) {
    hit_counts <- pubmed_hits[hits_up]
    hit_counts[is.na(hit_counts)] <- 0L
    ord    <- order(-hit_counts)[seq_len(min(10L, length(hits)))]
    hits_s <- hits[ord]
    ad_str <- paste0(hits_s, collapse = ", ")
  } else {
    ad_str <- "—"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    AD_genes     = ad_str,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(celltypes, build_row)) |>
  arrange(Celltype) |>
  dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

# ── Format with gt ─────────────────────────────────────────────────────────────
GT_SCALE <- 0.60  # uniform scale factor for PDF table size (1 = original)

gt_table <- summary_df |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 3 GO gene sets (by p-value)",
    AD_genes     = "AD-associated genes"
  ) |>
  fmt_markdown(columns = AD_genes) |>
  fmt(columns = Top_genesets, fns = function(x) x) |>
  cols_width(
    Celltype     ~ px(60  * GT_SCALE),
    N_genes      ~ px(45  * GT_SCALE),
    Top_genesets ~ px(450 * GT_SCALE),
    AD_genes     ~ px(262 * GT_SCALE)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_text(align = "center"),
    locations = list(cells_body(columns = N_genes),
                     cells_column_labels(columns = N_genes))
  ) |>
  tab_options(
    table.font.names      = c("Arial", "Helvetica", "sans-serif"),
    table.font.size       = 11 * GT_SCALE,
    data_row.padding      = px(6 * GT_SCALE),
    column_labels.padding = px(8 * GT_SCALE)
  )

out_path <- file.path(save_dir, "panel_G_gsea_summary_table.html")
gtsave(gt_table, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table_dfc.pdf")
tryCatch({
  gtsave(gt_table, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
})

bib_md   <- file.path(save_dir, "gsea_summary_bibliography_dfc.md")
bib_docx <- file.path(save_dir, "gsea_summary_bibliography_dfc.docx")
tryCatch({
  rmarkdown::pandoc_convert(bib_md, to = "docx", output = bib_docx)
  message("Saved: ", bib_docx)
}, error = function(e) {
  message("Bibliography docx conversion failed (requires pandoc): ", conditionMessage(e))
})

# ── FDR-corrected version (DFC) ─────────────────────────────────────────────
build_row_fdr <- function(ct) {
  genes   <- g_overlaps$Gene[g_overlaps$Celltype == ct]
  n_genes <- length(genes)

  ct_idx  <- which(celltypes == ct)
  pvals   <- gsea_out_fdr[[3 + ct_idx]]
  go_mask <- grepl("^GO", gsea_out_fdr$SetName)
  if (!is.null(pvals)) {
    go_pvals <- ifelse(go_mask, pvals, NA)
    sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
    n_sig    <- sum(sig_mask)
    if (n_sig == 0) {
      top_str <- "—"
    } else {
      top_idx  <- order(go_pvals)[1:min(3, n_sig)]
      top_sets <- gsea_out_fdr$SetName[top_idx]
      top_pval <- pvals[top_idx]
      top_str  <- paste0(
        seq_along(top_sets), ". ", top_sets,
        " (padj=", formatC(top_pval, format = "e", digits = 2), ")",
        collapse = "<br>"
      )
    }
  }

  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(ad_db_dfc)]
  hits_up     <- genes_upper[genes_upper %in% names(ad_db_dfc)]

  if (length(hits) > 0) {
    hit_counts <- pubmed_hits[hits_up]
    hit_counts[is.na(hit_counts)] <- 0L
    ord    <- order(-hit_counts)[seq_len(min(10L, length(hits)))]
    hits_s <- hits[ord]
    ad_str <- paste0(hits_s, collapse = ", ")
  } else {
    ad_str <- "—"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    AD_genes     = ad_str,
    stringsAsFactors = FALSE
  )
}

summary_df_fdr <- do.call(rbind, lapply(celltypes, build_row_fdr)) |>
  arrange(Celltype) |>
  dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

gt_table_fdr <- summary_df_fdr |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 3 GO gene sets (by FDR-adjusted p-value)",
    AD_genes     = "AD-associated genes"
  ) |>
  fmt_markdown(columns = AD_genes) |>
  fmt(columns = Top_genesets, fns = function(x) x) |>
  cols_width(
    Celltype     ~ px(60  * GT_SCALE),
    N_genes      ~ px(45  * GT_SCALE),
    Top_genesets ~ px(450 * GT_SCALE),
    AD_genes     ~ px(262 * GT_SCALE)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_text(align = "center"),
    locations = list(cells_body(columns = N_genes),
                     cells_column_labels(columns = N_genes))
  ) |>
  tab_options(
    table.font.names      = c("Arial", "Helvetica", "sans-serif"),
    table.font.size       = 11 * GT_SCALE,
    data_row.padding      = px(6 * GT_SCALE),
    column_labels.padding = px(8 * GT_SCALE)
  )

out_path <- file.path(save_dir, "panel_G_gsea_summary_table_dfc_fdr.html")
gtsave(gt_table_fdr, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table_dfc_fdr.pdf")
tryCatch({
  gtsave(gt_table_fdr, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
})


##########
# do mtg
g_overlaps <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.2/mtg_overlaps.csv"))

# Run GSEA
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func_optimized.R"))
b_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = T)
g_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = F)
gsea_out <- rbind(b_out, g_out)

# Calculate FDR
gsea_pval_to_fdr <- function(gsea_out) {                              
    meta_cols <- c("SetID", "SetName", "SetSize")                       
    pval_cols <- setdiff(colnames(gsea_out), meta_cols)                 
    fdr_out <- gsea_out                                                 
    fdr_out[, pval_cols] <- lapply(gsea_out[, pval_cols, drop = FALSE], 
                                    p.adjust, method = "BH")            
    fdr_out                                                             
  }                                                                     
gsea_out_fdr <- gsea_pval_to_fdr(gsea_out)

# Significant genesets per celltype (p < GSEA_PTHRESH)
# GSEA_PTHRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   gsea_out$SetName[!is.na(pvals) & pvals < GSEA_PTHRESH]
# }), celltypes)

# ── AD gene database ────────────────────────────────────────────────────────
library(gt)

mtg_hits    <- fread(file.path(save_dir, "ad_db_summary_table_mtg.csv"))
pubmed_hits <- setNames(as.integer(mtg_hits$Pubmed_total_hits), mtg_hits$Gene)

# ── Build table ─────────────────────────────────────────────────
celltypes  <- unique(g_overlaps$Celltype)

# Significant genesets per celltype (BH-adjusted p < GSEA_FDR_THRESH)
# GSEA_FDR_THRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   fdr   <- p.adjust(pvals, method = "BH")
#   gsea_out$SetName[!is.na(fdr) & fdr < GSEA_FDR_THRESH]
# }), celltypes)
# qsave(gsea_sig, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/gsea_sig_mtg.qs"))

build_row <- function(ct) {
  genes  <- g_overlaps$Gene[g_overlaps$Celltype == ct]
  n_genes <- length(genes)

  # Top 3 GO gene sets for this celltype (columns assumed in same order as celltypes)
  ct_idx  <- which(celltypes == ct)
  pvals   <- gsea_out[[3 + ct_idx]]
  go_mask <- grepl("^GO", gsea_out$SetName)
  if (!is.null(pvals)) {
    go_pvals <- ifelse(go_mask, pvals, NA)
    sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
    n_sig    <- sum(sig_mask)
    if (n_sig == 0) {
      top_str <- "—"
    } else {
      top_idx  <- order(go_pvals)[1:min(3, n_sig)]
      top_sets <- gsea_out$SetName[top_idx]
      top_pval <- pvals[top_idx]
      top_str  <- paste0(
        seq_along(top_sets), ". ", top_sets,
        " (p=", formatC(top_pval, format = "e", digits = 2), ")",
        collapse = "<br>"
      )
    }
  }

  # AD genes: top 10 by PubMed hits
  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(ad_db_mtg)]
  hits_up     <- genes_upper[genes_upper %in% names(ad_db_mtg)]

  if (length(hits) > 0) {
    hit_counts <- pubmed_hits[hits_up]
    hit_counts[is.na(hit_counts)] <- 0L
    ord    <- order(-hit_counts)[seq_len(min(10L, length(hits)))]
    hits_s <- hits[ord]
    ad_str <- paste0(hits_s, collapse = ", ")
  } else {
    ad_str <- "—"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    AD_genes     = ad_str,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(celltypes, build_row)) |>
  arrange(Celltype) |>
  dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

# ── Format with gt ─────────────────────────────────────────────────────────────
GT_SCALE <- 0.60  # uniform scale factor for PDF table size (1 = original)

gt_table <- summary_df |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 3 GO gene sets (by p-value)",
    AD_genes     = "AD-associated genes"
  ) |>
  fmt_markdown(columns = AD_genes) |>
  fmt(columns = Top_genesets, fns = function(x) x) |>
  cols_width(
    Celltype     ~ px(60  * GT_SCALE),
    N_genes      ~ px(45  * GT_SCALE),
    Top_genesets ~ px(450 * GT_SCALE),
    AD_genes     ~ px(262 * GT_SCALE)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_text(align = "center"),
    locations = list(cells_body(columns = N_genes),
                     cells_column_labels(columns = N_genes))
  ) |>
  tab_options(
    table.font.names      = c("Arial", "Helvetica", "sans-serif"),
    table.font.size       = 11 * GT_SCALE,
    data_row.padding      = px(6 * GT_SCALE),
    column_labels.padding = px(8 * GT_SCALE)
  )

out_path <- file.path(save_dir, "panel_G_gsea_summary_table_mtg.html")
gtsave(gt_table, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table_mtg.pdf")
tryCatch({
  gtsave(gt_table, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
})

bib_md   <- file.path(save_dir, "gsea_summary_bibliography_mtg.md")
bib_docx <- file.path(save_dir, "gsea_summary_bibliography_mtg.docx")
tryCatch({
  rmarkdown::pandoc_convert(bib_md, to = "docx", output = bib_docx)
  message("Saved: ", bib_docx)
}, error = function(e) {
  message("Bibliography docx conversion failed (requires pandoc): ", conditionMessage(e))
})

# ── FDR-corrected version (MTG) ─────────────────────────────────────────────
build_row_fdr <- function(ct) {
  genes   <- g_overlaps$Gene[g_overlaps$Celltype == ct]
  n_genes <- length(genes)

  ct_idx  <- which(celltypes == ct)
  pvals   <- gsea_out_fdr[[3 + ct_idx]]
  go_mask <- grepl("^GO", gsea_out_fdr$SetName)
  if (!is.null(pvals)) {
    go_pvals <- ifelse(go_mask, pvals, NA)
    sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
    n_sig    <- sum(sig_mask)
    if (n_sig == 0) {
      top_str <- "—"
    } else {
      top_idx  <- order(go_pvals)[1:min(3, n_sig)]
      top_sets <- gsea_out_fdr$SetName[top_idx]
      top_pval <- pvals[top_idx]
      top_str  <- paste0(
        seq_along(top_sets), ". ", top_sets,
        " (padj=", formatC(top_pval, format = "e", digits = 2), ")",
        collapse = "<br>"
      )
    }
  }

  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(ad_db_mtg)]
  hits_up     <- genes_upper[genes_upper %in% names(ad_db_mtg)]

  if (length(hits) > 0) {
    hit_counts <- pubmed_hits[hits_up]
    hit_counts[is.na(hit_counts)] <- 0L
    ord    <- order(-hit_counts)[seq_len(min(10L, length(hits)))]
    hits_s <- hits[ord]
    ad_str <- paste0(hits_s, collapse = ", ")
  } else {
    ad_str <- "—"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    AD_genes     = ad_str,
    stringsAsFactors = FALSE
  )
}

summary_df_fdr <- do.call(rbind, lapply(celltypes, build_row_fdr)) |>
  arrange(Celltype) |>
  dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

gt_table_fdr <- summary_df_fdr |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 3 GO gene sets (by FDR-adjusted p-value)",
    AD_genes     = "AD-associated genes"
  ) |>
  fmt_markdown(columns = AD_genes) |>
  fmt(columns = Top_genesets, fns = function(x) x) |>
  cols_width(
    Celltype     ~ px(60  * GT_SCALE),
    N_genes      ~ px(45  * GT_SCALE),
    Top_genesets ~ px(450 * GT_SCALE),
    AD_genes     ~ px(262 * GT_SCALE)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_text(align = "center"),
    locations = list(cells_body(columns = N_genes),
                     cells_column_labels(columns = N_genes))
  ) |>
  tab_options(
    table.font.names      = c("Arial", "Helvetica", "sans-serif"),
    table.font.size       = 11 * GT_SCALE,
    data_row.padding      = px(6 * GT_SCALE),
    column_labels.padding = px(8 * GT_SCALE)
  )

out_path <- file.path(save_dir, "panel_G_gsea_summary_table_mtg_fdr.html")
gtsave(gt_table_fdr, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table_mtg_fdr.pdf")
tryCatch({
  gtsave(gt_table_fdr, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
})

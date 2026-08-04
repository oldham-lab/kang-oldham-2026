library(tidyverse)
library(qs)
library(data.table)
library(openxlsx)

# Table S4: Gene expression modeling results for all individuals, regions, and technology platforms
# (Mean Expression Percentile (MEP), Adj. R2, and RMSE).
#
# Output: Two .xlsx files — one for Jorstad et al., one for SEAAD 2024 (Gabitto et al.)
# Each file has one tab per unique combination of Donor x Platform x Region x Celltype,
# plus a Legend tab summarising columns and the combinations present.

# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s6_7/v2/"))
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)


# ============================================================
# Helper: calculate percentile
# ============================================================
calc_per <- function(x) {
  (1 - (rank(-x, na.last = "keep") / max(rank(-x, na.last = "keep")))) * 100
}


# ============================================================
# Table creation function
# ============================================================
create_table <- function(homedir, cell_anno, genemap, donornames,
                         san_mean, platform, region) {

  plotdflist_sub   <- qread(paste0(homedir, "/fig1bdflist_subclass.qs"))
  plotdflist_super <- qread(paste0(homedir, "/fig1bdflist_supertype.qs"))

  san_mean <- lapply(san_mean, function(x) {
    as.data.frame(x) |>
      rename(mean = "x") |>
      rownames_to_column(var = "Gene")
  })

  san_mean_per <- lapply(san_mean, \(x) {
    x |>
      mutate(mean_expr_pcntile = calc_per(mean)) |>
      select(Gene, mean_expr_pcntile)
  })

  # Collect gene symbols from each donor (pcnt.var = 0)
  genesymlist <- lapply(seq_along(donornames), function(i) {
    exdir     <- list.files(paste0(homedir, "/donor", i, "/SyntheticDatasets/"), full.names = TRUE)
    exdirexpr <- exdir[grep("EXPRLIST", exdir)]
    exdirexpr <- exdirexpr[grep("0pcntVar", exdirexpr)]
    obj <- readRDS(exdirexpr)
    obj[[1]][, 2]
  })

  # Add gene symbols and filter to pcnt.var == 0
  plotdflist_sub <- mapply(function(x, y) {
    x |> dplyr::filter(pcnt.var == 0) |> mutate(Gene = y)
  }, plotdflist_sub, genesymlist, SIMPLIFY = FALSE)

  plotdflist_super <- mapply(function(x, y) {
    x |> dplyr::filter(pcnt.var == 0) |> mutate(Gene = y)
  }, plotdflist_super, genesymlist, SIMPLIFY = FALSE)

  # Replace mean with SN-derived mean expression percentile; add metadata columns
  plotdflist_sub <- mapply(function(a, b, c) {
    a |>
      dplyr::filter(pcnt.var == 0) |>
      dplyr::select(!mean) |>
      left_join(b, by = "Gene") |>
      mutate(Donor = c, Platform = platform, Region = region, Celltype = "subclass") |>
      select(Gene, Donor, Platform, Region, Celltype, mean_expr_pcntile, adj_r2, rmse) |>
      arrange(Gene)
  }, plotdflist_sub, san_mean_per, donornames, SIMPLIFY = FALSE) |>
    do.call(what = "rbind")

  plotdflist_super <- mapply(function(a, b, c) {
    a |>
      dplyr::filter(pcnt.var == 0) |>
      dplyr::select(!mean) |>
      left_join(b, by = "Gene") |>
      mutate(Donor = c, Platform = platform, Region = region, Celltype = "supertype") |>
      select(Gene, Donor, Platform, Region, Celltype, mean_expr_pcntile, adj_r2, rmse) |>
      arrange(Gene)
  }, plotdflist_super, san_mean_per, donornames, SIMPLIFY = FALSE) |>
    do.call(what = "rbind")

  rbind(plotdflist_sub, plotdflist_super)
}


# ============================================================
# Data generation — Jorstad et al.
# ============================================================

# DFC — Chromium v3
homedir    <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/")
cell_anno  <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
genemap    <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"), data.table = FALSE)
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
san_mean   <- qread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data_other/lein_dfc/gene_count_means_byDonor.qs"))

jor_dfc <- create_table(homedir, cell_anno, genemap, donornames, san_mean,
                        platform = "Cv3", region = "DLPFC")

# MTG — Chromium v3
homedir    <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/MTG_indiv_donor/")
cell_anno  <- fread(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_MTG.csv"), data.table = FALSE)
genemap    <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"), data.table = FALSE)
donornames <- c("H200.1023", "H200.1025", "H200.1030")
san_mean   <- qread(file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_10x/gene_count_means_byDonor.qs"))

jor_mtg <- create_table(homedir, cell_anno, genemap, donornames, san_mean,
                        platform = "Cv3", region = "MTG")

# V1 — Chromium v3
homedir    <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/V1_indiv_donor/")
cell_anno  <- fread(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_V1.csv"), data.table = FALSE)
cell_anno  <- cell_anno[, c(3, 1, 2, 4:8)]
genemap    <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"), data.table = FALSE)
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
san_mean   <- qread(file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_v1_10x/gene_count_means_byDonor.qs"))

jor_v1 <- create_table(homedir, cell_anno, genemap, donornames, san_mean,
                        platform = "Cv3", region = "V1")

# DLPFC — SMART-seq v4
homedir    <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_SSv4_indiv_donor/")
cell_anno  <- fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_ssv4/lein_dfc_metadata.csv"), data.table = FALSE)
genemap    <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"), data.table = FALSE)
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
san_mean   <- qread(file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_ssv4/gene_count_means_byDonor.qs"))

jor_ssv4_dfc <- create_table(homedir, cell_anno, genemap, donornames, san_mean,
                              platform = "SSv4", region = "DLPFC")

# MTG — SMART-seq v4
homedir    <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_SSv4_indiv_donor_mtg/")
cell_anno  <- fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_ssv4/lein_mtg_metadata.csv"), data.table = FALSE)
genemap    <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"), data.table = FALSE)
donornames <- c("H200.1023", "H200.1025", "H200.1030")
san_mean   <- qread(file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_ssv4/gene_count_means_byDonor.qs"))

jor_ssv4_mtg <- create_table(homedir, cell_anno, genemap, donornames, san_mean,
                              platform = "SSv4", region = "MTG")

# Combine all Jorstad results
jor_out <- rbind(jor_mtg, jor_dfc, jor_v1, jor_ssv4_mtg, jor_ssv4_dfc)


# ============================================================
# Data generation — SEAAD 2024 (Gabitto et al.)
# ============================================================
homedir   <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/SEAAD2024con_indiv_donor")
cell_anno <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv"), data.table = FALSE) |>
  dplyr::filter(`Overall AD neuropathological Change` == "Not AD")
donornames <- unique(cell_anno$`Donor ID`)
san_mean   <- qread(file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_con/gene_count_means_byDonor.qs"))

sea_dfc <- create_table(homedir, cell_anno, genemap, donornames, san_mean,
                        platform = "Cv3", region = "DLPFC")


# ============================================================
# openxlsx formatting helpers
# ============================================================
SIG_FIGS <- 4

round_sig_df <- function(df, n = SIG_FIGS) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) signif(col, n) else col
  })
  df
}

# Truncate tab name to Excel's 31-character limit
make_tab_name <- function(donor, platform, region, celltype) {
  nm <- paste(donor, platform, region, celltype, sep = " | ")
  substr(nm, 1, 31)
}

format_data_sheet <- function(wb, sheet, df) {
  n_rows <- nrow(df)
  n_cols <- ncol(df)

  header_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    halign         = "center",
    valign         = "center",
    fgFill         = "#D9E1F2",
    border         = "Bottom",
    borderColour   = "#4472C4",
    borderStyle    = "medium"
  )
  body_style <- createStyle(
    fontName = "Arial",
    fontSize = 11,
    halign   = "center",
    valign   = "center"
  )

  addStyle(wb, sheet, header_style,
           rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  if (n_rows > 0)
    addStyle(wb, sheet, body_style,
             rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)

  # Size each column to the widest of its header / values (+ padding) so nothing clips
  col_widths <- vapply(seq_len(n_cols), function(j) {
    max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm = TRUE)
  }, numeric(1))
  setColWidths(wb, sheet, cols = 1:n_cols, widths = col_widths + 2)
  freezePane(wb, sheet, firstRow = TRUE)
}

write_legend_section <- function(wb, sheet, title, legend_df, start_row) {
  n_cols <- ncol(legend_df)

  # Blue section header spanning all columns of legend_df
  section_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    fontColour     = "#FFFFFF",
    fgFill         = "#4472C4",
    halign         = "left",
    valign         = "center"
  )
  mergeCells(wb, sheet, cols = 1:n_cols, rows = start_row)
  writeData(wb, sheet, x = title, startRow = start_row, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet, section_style, rows = start_row, cols = 1:n_cols, gridExpand = TRUE)

  # Sub-header row spanning all columns
  subheader_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    fgFill         = "#D9E1F2",
    halign         = "left",
    valign         = "center",
    border         = "Bottom",
    borderColour   = "#4472C4",
    borderStyle    = "medium"
  )
  writeData(wb, sheet, x = legend_df, startRow = start_row + 1, startCol = 1, colNames = TRUE)
  addStyle(wb, sheet, subheader_style,
           rows = start_row + 1, cols = 1:n_cols, gridExpand = TRUE)

  # Body rows spanning all columns
  body_style <- createStyle(
    fontName = "Arial",
    fontSize = 11,
    halign   = "left",
    valign   = "center",
    wrapText = TRUE
  )
  body_rows <- (start_row + 2):(start_row + 1 + nrow(legend_df))
  addStyle(wb, sheet, body_style, rows = body_rows, cols = 1:n_cols, gridExpand = TRUE)

  # Return next available row (leave one blank row gap)
  start_row + 1 + nrow(legend_df) + 2
}


# ============================================================
# Core workbook builder
# ============================================================
build_workbook <- function(data, filename, dataset_label, title) {

  wb <- createWorkbook()
  modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

  # Split into one data frame per unique Donor x Platform x Region x Celltype combination
  combos <- data |>
    distinct(Donor, Platform, Region, Celltype) |>
    arrange(Region, Platform, Donor, Celltype)

  for (k in seq_len(nrow(combos))) {
    combo    <- combos[k, ]
    tab_df   <- data |>
      dplyr::filter(Donor    == combo$Donor,
                    Platform == combo$Platform,
                    Region   == combo$Region,
                    Celltype == combo$Celltype) |>
      round_sig_df()

    tab_name <- make_tab_name(combo$Donor, combo$Platform, combo$Region, combo$Celltype)
    addWorksheet(wb, tab_name)
    writeData(wb, tab_name, tab_df, rowNames = FALSE)
    format_data_sheet(wb, tab_name, tab_df)
  }

  # ---- Legend tab ----
  addWorksheet(wb, "Legend")

  TITLE_OFFSET <- 2  # title row + blank spacer row before legend content

  # Title row
  writeData(wb, "Legend", data.frame(x = title),
            rowNames = FALSE, colNames = FALSE, startRow = 1)
  mergeCells(wb, "Legend", cols = 1:2, rows = 1)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold",
                       halign = "left", valign = "center", wrapText = TRUE),
           rows = 1, cols = 1:2, gridExpand = TRUE)

  col_legend <- data.frame(
    Column = c("Gene", "Donor", "Platform", "Region", "Celltype",
               "mean_expr_pcntile", "adj_r2", "rmse"),
    Description = c(
      "HGNC gene symbol",
      "Donor ID",
      "Sequencing platform used to generate the data (Cv3 = 10x Chromium v3; SSv4 = SMART-seq v4)",
      "Brain region from which nuclei were collected",
      "Level of cell type annotation used in the model (subclass or supertype)",
      "Mean expression percentile of the gene across all cells; higher values indicate more highly expressed genes",
      "Adjusted R-squared from the gene expression model; reflects the proportion of variance in expression explained by cell type composition",
      "Root mean square error from the gene expression model; reflects the average prediction error in expression units"
    ),
    stringsAsFactors = FALSE
  )

  # Build overview table of all combinations present
  overview <- data |>
    distinct(Donor, Platform, Region, Celltype) |>
    arrange(Region, Platform, Donor, Celltype)

  next_row <- write_legend_section(wb, "Legend",
                                   title      = "Column descriptions",
                                   legend_df  = col_legend,
                                   start_row  = 1 + TITLE_OFFSET)

  next_row <- write_legend_section(wb, "Legend",
                                   title      = paste0("Combinations present — ", dataset_label),
                                   legend_df  = overview,
                                   start_row  = next_row)

  # Note on tab naming
  note_style <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "italic")
  writeData(wb, "Legend",
            x         = "Note: Each tab is labelled as Donor | Platform | Region | Celltype. Tab names are truncated to 31 characters (Excel limit).",
            startRow  = next_row,
            startCol  = 1,
            colNames  = FALSE)
  addStyle(wb, "Legend", note_style, rows = next_row, cols = 1:2, gridExpand = TRUE)

  # Format Legend sheet
  setColWidths(wb, "Legend", cols = 1, widths = 25)
  setColWidths(wb, "Legend", cols = 2, widths = 80)
  freezePane(wb, "Legend", firstRow = TRUE)

  # Legend tab first
  leg_idx <- which(names(wb) == "Legend")
  worksheetOrder(wb) <- c(leg_idx, setdiff(seq_along(names(wb)), leg_idx))

  out_path <- file.path(save_dir, filename)
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message("Saved: ", out_path)
}


# ============================================================
# Build and save workbooks
# ============================================================
build_workbook(jor_out, filename = "s6_jorstad.xlsx",  dataset_label = "Jorstad et al.",
               title = "Table S6 | Pseudobulk modeling of gene expression as a function of cell-type abundance (Jorstad et al.)")
build_workbook(sea_dfc, filename = "s7_seaad2024.xlsx", dataset_label = "SEAAD 2024 (Gabitto et al.)",
               title = "Table S7: Pseudobulk modeling of gene expression as a function of cell-type abundance (Gabitto et al.)")
